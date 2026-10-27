#!/usr/bin/env Rscript


options(stringsAsFactors = FALSE, warn = 1)
set.seed(20260922)

# ------------------------------ paths -----------------------------------------
get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  hit <- grep("^--file=", args, value = TRUE)
  if (length(hit) == 0) return(normalizePath(getwd(), mustWork = TRUE))
  normalizePath(dirname(sub("^--file=", "", hit[[1]])), mustWork = TRUE)
}
script_dir <- get_script_dir()
ROOT <- if (basename(script_dir) == "supplementary") dirname(script_dir) else script_dir
setwd(ROOT)


FILES <- list(
  per_genome = "per_genome.csv",
  per_inversion = "per_inversion.csv",
  animal = "phylogenetic/animal_mcmctree_augmented_verified.nwk",
  plant =  "phylogenetic/plant_mcmctree_sanitized.nwk")

OUTDIR <- file.path(ROOT, "phylogenetic_revision_FINAL")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

# ------------------------------ packages --------------------------------------

suppressPackageStartupMessages({
  library(ape)
  library(phytools)
  library(phylolm)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(ggrepel)
})

# ------------------------------ style -----------------------------------------
KINGDOM_COLS <- c(Metazoa = "#E1812C", Viridiplantae = "#00A087")
TRAIT_LABELS <- c(
  log10_inv_count_per_Gb = "Inversion count per Gb",
  inv_count_total="Total inversion count",
  median_log10_inv_length = "Median inversion length",
  logit_prop_large = "Proportion >1 Mb",
  log10_median_pi_ratio = "Median pi_inv / pi_genome",
  gene_density = "Gene density"
)
PRED_LABELS <- c(
  z_log10_genome_size = "Genome size",
  z_TE_total_pct = "TE content",
  z_pi_genome_wide = "Genome-wide pi",
  z_gene_density = "Gene density"
)

theme_nature <- function(base_size = 8) {
  theme_classic(base_size = base_size) +
    theme(
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black"),
      axis.line = element_line(linewidth = 0.35, color = "black"),
      axis.ticks = element_line(linewidth = 0.35, color = "black"),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", color = "black"),
      legend.title = element_blank(),
      legend.key.height = grid::unit(0.35, "cm"),
      plot.tag = element_text(face = "bold", size = base_size + 1),
      plot.margin = margin(4, 5, 4, 5)
    )
}

# ------------------------------ helpers ---------------------------------------
as_num <- function(x) suppressWarnings(as.numeric(x))

safe_median <- function(x) {
  x <- as_num(x)
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  median(x)
}

safe_fraction_gt <- function(x, cutoff) {
  x <- as_num(x)
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  mean(x > cutoff)
}

clamp01 <- function(x, eps = 1e-6) pmin(pmax(x, eps), 1 - eps)

zscore <- function(x) {
  x <- as_num(x)
  s <- sd(x, na.rm = TRUE)
  m <- mean(x, na.rm = TRUE)
  if (!is.finite(s) || s <= 0) stop("Cannot z-standardize a constant/non-finite variable.", call. = FALSE)
  (x - m) / s
}

require_columns <- function(df, cols, label) {
  miss <- setdiff(cols, names(df))
  if (length(miss) > 0) {
    stop(label, " is missing required columns: ", paste(miss, collapse = ", "), call. = FALSE)
  }
}

read_tree_strict <- function(path, label) {
  tr <- ape::read.tree(file = path)
  if (!inherits(tr, "phylo")) stop("Cannot read tree: ", label, call. = FALSE)
  tr$tip.label <- trimws(tr$tip.label)
  if (length(tr$tip.label) < 4) stop("Tree has <4 tips: ", label, call. = FALSE)
  if (anyNA(tr$tip.label) || any(tr$tip.label == "")) stop("Empty/NA tip label: ", label, call. = FALSE)
  if (anyDuplicated(tr$tip.label)) stop("Duplicated tip label: ", label, call. = FALSE)
  if (is.null(tr$edge.length)) stop("Missing branch lengths: ", label, call. = FALSE)
  if (any(!is.finite(tr$edge.length))) stop("Non-finite branch length: ", label, call. = FALSE)
  if (any(tr$edge.length < -1e-10)) stop("Negative branch length: ", label, call. = FALSE)

  # Tiny numerical deviations from ultrametricity are expected in MCMCTree output.
  dep <- ape::node.depth.edgelength(tr)[seq_len(ape::Ntip(tr))]
  rel_spread <- (max(dep) - min(dep)) / max(dep)
  if (rel_spread > 1e-4) {
    warning(label, " is not closely ultrametric (relative root-to-tip spread = ",
            signif(rel_spread, 4), ").")
  }
  tr
}

tree_audit <- function(tr, label) {
  dep <- ape::node.depth.edgelength(tr)[seq_len(ape::Ntip(tr))]
  tibble(
    tree = label,
    n_tips = ape::Ntip(tr),
    n_nodes = ape::Nnode(tr),
    duplicate_tips = sum(duplicated(tr$tip.label)),
    negative_edges = sum(tr$edge.length < -1e-10),
    zero_edges = sum(abs(tr$edge.length) <= 1e-12),
    root_to_tip_min = min(dep),
    root_to_tip_max = max(dep),
    root_to_tip_spread = max(dep) - min(dep),
    relative_spread = (max(dep) - min(dep)) / max(dep)
  )
}

prune_and_order <- function(tree, dat) {
  if (!"tip" %in% names(dat)) stop("Analysis data must contain column 'tip'.", call. = FALSE)
  if (anyDuplicated(dat$tip)) stop("Duplicate tip labels in analysis data.", call. = FALSE)

  keep <- intersect(tree$tip.label, dat$tip)
  if (length(keep) < 4) stop("Fewer than four matched tips after pruning.", call. = FALSE)

  tr <- ape::keep.tip(tree, keep)
  d <- dat[match(tr$tip.label, dat$tip), , drop = FALSE]
  if (anyNA(d$tip)) stop("Tree/data matching failure.", call. = FALSE)
  if (!identical(tr$tip.label, d$tip)) stop("Tree/data ordering failure.", call. = FALSE)
  rownames(d) <- d$tip
  list(tree = tr, data = d)
}

map_tree_to_genomes <- function(tree, genome_df, kingdom_name) {
  dat <- genome_df %>% filter(.data$kingdom == kingdom_name)
  tips <- tree$tip.label

  # Priority 1: exact sample_id.
  exact_idx <- match(tips, dat$sample_id)
  exact_ok <- !is.na(exact_idx)
  exact <- tibble(
    tip = tips[exact_ok],
    row_id = dat$.row_id[exact_idx[exact_ok]],
    mapping_method = "sample_id_exact"
  )

  # Priority 2: exact unique species label after spaces -> underscores.
  used_rows <- exact$row_id
  remaining_tips <- setdiff(tips, exact$tip)

  sp_counts <- dat %>% count(.data$species_tip, name = "n_rows")
  unique_sp <- dat %>%
    left_join(sp_counts, by = "species_tip") %>%
    filter(.data$n_rows == 1, !(.data$.row_id %in% used_rows))

  idx <- match(remaining_tips, unique_sp$species_tip)
  ok <- !is.na(idx)
  fallback <- tibble(
    tip = remaining_tips[ok],
    row_id = unique_sp$.row_id[idx[ok]],
    mapping_method = "unique_species_fallback"
  )

  mp <- bind_rows(exact, fallback)
  if (anyDuplicated(mp$tip)) stop("A tree tip was mapped more than once.", call. = FALSE)
  if (anyDuplicated(mp$row_id)) stop("A genome row was mapped to >1 tree tip.", call. = FALSE)

  mapped <- mp %>% left_join(dat, by = c("row_id" = ".row_id"))
  unmapped_tree_tips <- setdiff(tips, mapped$tip)
  unmapped_genomes <- dat %>% filter(!(.data$.row_id %in% mapped$row_id))

  list(mapped = mapped,
       unmapped_tree_tips = unmapped_tree_tips,
       unmapped_genomes = unmapped_genomes)
}

# -------------------- custom phylogenetic variance partition ------------------
# Model: y ~ intercept + phylogenetic random covariance + independent residual
# V = h*C + (1-h)*I, where C is scaled to mean(diag(C))=1.
# h is therefore an ML estimate of the fraction of variance assigned to the
# phylogenetic covariance component under this fitted model.

profile_h2 <- function(tree, dat, trait) {
  use <- dat %>%
    select(tip, all_of(trait)) %>%
    filter(is.finite(.data[[trait]]))

  po <- prune_and_order(tree, use)
  y <- as.numeric(po$data[[trait]])
  n <- length(y)
  if (n < 8 || sd(y) <= 0) {
    return(tibble(trait = trait, n = n, H2_phy = NA_real_,
                  ci_low = NA_real_, ci_high = NA_real_,
                  logLik = NA_real_))
  }

  C <- ape::vcv.phylo(po$tree, corr = FALSE)
  C <- C[po$data$tip, po$data$tip, drop = FALSE]
  C <- C / mean(diag(C))
  I <- diag(n)
  one <- rep(1, n)

  nll <- function(h) {
    if (!is.finite(h) || h < 0 || h > 0.999999) return(Inf)
    V <- h * C + (1 - h) * I
    R <- try(chol(V), silent = TRUE)
    if (inherits(R, "try-error")) return(Inf)

    Vinv_y <- backsolve(R, forwardsolve(t(R), y))
    Vinv_1 <- backsolve(R, forwardsolve(t(R), one))
    beta <- sum(one * Vinv_y) / sum(one * Vinv_1)
    r <- y - beta
    Vinv_r <- backsolve(R, forwardsolve(t(R), r))
    rss <- sum(r * Vinv_r)
    sigma2 <- rss / n
    if (!is.finite(sigma2) || sigma2 <= 0) return(Inf)

    logdetV <- 2 * sum(log(diag(R)))
    0.5 * (n * log(2 * pi * sigma2) + logdetV + n)
  }

  opt <- optimize(nll, interval = c(0, 0.999999), tol = 1e-7)
  hhat <- opt$minimum
  min_nll <- opt$objective

  # 95% profile-likelihood interval: 2*DeltaLL <= chi^2_1(0.95)
  cutoff <- min_nll + 0.5 * qchisq(0.95, df = 1)
  grid <- unique(sort(c(seq(0, 0.999999, length.out = 401), hhat)))
  vals <- vapply(grid, nll, numeric(1))
  ok <- vals <= cutoff
  ci_low <- if (any(ok)) min(grid[ok]) else NA_real_
  ci_high <- if (any(ok)) max(grid[ok]) else NA_real_

  tibble(
    trait = trait, n = n,
    H2_phy = hhat,
    ci_low = ci_low,
    ci_high = ci_high,
    logLik = -min_nll
  )
}

# -------------------------- phylogenetic signal -------------------------------
run_phylosig_one <- function(tree, dat, trait, kingdom_name) {
  use <- dat %>%
    select(tip, all_of(trait)) %>%
    filter(is.finite(.data[[trait]]))

  po <- prune_and_order(tree, use)
  x <- po$data[[trait]]
  names(x) <- po$data$tip

  if (length(x) < 8 || sd(x) <= 0) {
    return(tibble(
      kingdom = kingdom_name, trait = trait, n = length(x),
      lambda = NA_real_, lambda_P = NA_real_,
      K = NA_real_, K_P = NA_real_
    ))
  }

  lam <- phytools::phylosig(po$tree, x, method = "lambda", test = TRUE)
  kval <- phytools::phylosig(po$tree, x, method = "K", test = TRUE, nsim = 999)

  tibble(
    kingdom = kingdom_name,
    trait = trait,
    n = length(x),
    lambda = as.numeric(lam$lambda),
    lambda_P = as.numeric(lam$P),
    K = as.numeric(kval$K),
    K_P = as.numeric(kval$P)
  )
}

# ------------------------------- PGLS -----------------------------------------
run_pgls <- function(tree, dat, response, kingdom_name,
                     include_gene_density = FALSE,
                     model_label = "primary") {

  pred_raw <- c("log10_genome_size", "TE_total_pct", "pi_genome_wide")
  if (include_gene_density) pred_raw <- c(pred_raw, "gene_density")

  needed <- c("tip", response, pred_raw)
  use <- dat %>%
    select(all_of(needed)) %>%
    filter(if_all(all_of(c(response, pred_raw)), ~ !is.na(.x))) %>%
    filter(if_all(all_of(c(response, pred_raw)), ~ is.finite(.x)))

  po <- prune_and_order(tree, use)
  d <- po$data

  d$z_response <- zscore(d[[response]])
  d$z_log10_genome_size <- zscore(d$log10_genome_size)
  d$z_TE_total_pct <- zscore(d$TE_total_pct)
  d$z_pi_genome_wide <- zscore(d$pi_genome_wide)
  if (include_gene_density) d$z_gene_density <- zscore(d$gene_density)
  rownames(d) <- d$tip

  rhs <- c("z_log10_genome_size", "z_TE_total_pct", "z_pi_genome_wide")
  if (include_gene_density) rhs <- c(rhs, "z_gene_density")
  fml <- as.formula(paste("z_response ~", paste(rhs, collapse = " + ")))

  fit <- phylolm::phylolm(
    formula = fml,
    data = d,
    phy = po$tree,
    model = "lambda",
    REML = FALSE
  )

  beta <- fit$coefficients
  se <- sqrt(diag(fit$vcov))
  z <- beta / se
  p <- 2 * pnorm(abs(z), lower.tail = FALSE)

  coef_tbl <- tibble(
    kingdom = kingdom_name,
    response = response,
    model = model_label,
    n = fit$n,
    term = names(beta),
    estimate = as.numeric(beta),
    std_error = as.numeric(se),
    ci_low = as.numeric(beta - 1.96 * se),
    ci_high = as.numeric(beta + 1.96 * se),
    P = as.numeric(p)
  )

  model_tbl <- tibble(
    kingdom = kingdom_name,
    response = response,
    model = model_label,
    n = fit$n,
    lambda = as.numeric(fit$optpar),
    logLik = as.numeric(fit$logLik),
    AIC = as.numeric(fit$aic),
    R2 = as.numeric(fit$r.squared),
    adjusted_R2 = as.numeric(fit$adj.r.squared)
  )

  list(coefficients = coef_tbl, model = model_tbl, data = d)
}

run_simple_pgls <- function(tree, dat, response, predictor, kingdom_name) {
  use <- dat %>%
    select(tip, all_of(c(response, predictor))) %>%
    filter(if_all(all_of(c(response, predictor)), ~ !is.na(.x))) %>%
    filter(if_all(all_of(c(response, predictor)), ~ is.finite(.x)))

  po <- prune_and_order(tree, use)
  d <- po$data
  d$z_response <- zscore(d[[response]])
  d$z_predictor <- zscore(d[[predictor]])
  rownames(d) <- d$tip

  fit <- phylolm::phylolm(
    z_response ~ z_predictor,
    data = d,
    phy = po$tree,
    model = "lambda",
    REML = FALSE
  )

  b <- fit$coefficients["z_predictor"]
  se <- sqrt(diag(fit$vcov))["z_predictor"]
  p <- 2 * pnorm(abs(b / se), lower.tail = FALSE)

  tibble(
    kingdom = kingdom_name,
    response = response,
    predictor = predictor,
    n = fit$n,
    beta_standardized = as.numeric(b),
    std_error = as.numeric(se),
    ci_low = as.numeric(b - 1.96 * se),
    ci_high = as.numeric(b + 1.96 * se),
    P = as.numeric(p),
    lambda = as.numeric(fit$optpar),
    R2 = as.numeric(fit$r.squared)
  )
}

# ------------------------------ VIF --------------------------------------------
calc_vif <- function(df, predictors, kingdom_name, response, model_label) {
  d <- df %>%
    select(all_of(predictors)) %>%
    filter(if_all(everything(), ~ !is.na(.x))) %>%
    filter(if_all(everything(), ~ is.finite(.x)))

  out <- lapply(predictors, function(v) {
    others <- setdiff(predictors, v)
    if (length(others) == 0) return(NA_real_)
    fit <- lm(reformulate(others, response = v), data = d)
    r2 <- summary(fit)$r.squared
    if (!is.finite(r2) || r2 >= 1) Inf else 1 / (1 - r2)
  })

  tibble(
    kingdom = kingdom_name,
    response = response,
    model = model_label,
    predictor = predictors,
    VIF = unlist(out)
  )
}

# -------------------------- cross-kingdom PGLS --------------------------------
run_cross_kingdom_pgls <- function(tree, animal_dat, plant_dat,
                                   response, root_age_Ma) {

  use_cols <- c("tip", "kingdom", response,
                "log10_genome_size", "TE_total_pct", "pi_genome_wide","gene_density")

  all_dat <- bind_rows(
    animal_dat %>% select(all_of(use_cols)),
    plant_dat %>% select(all_of(use_cols))
  ) %>%
    filter(
      is.finite(.data[[response]]),
      is.finite(.data$log10_genome_size),
      is.finite(.data$TE_total_pct),
      is.finite(.data$pi_genome_wide),
      is.finite(.data$gene_density)
    )

  po <- prune_and_order(tree, all_dat)
  d <- po$data
  d$kingdom <- factor(d$kingdom, levels = c("Metazoa", "Viridiplantae"))
  d$z_response <- zscore(d[[response]])
  d$z_log10_genome_size <- zscore(d$log10_genome_size)
  d$z_TE_total_pct <- zscore(d$TE_total_pct)
  d$z_pi_genome_wide <- zscore(d$pi_genome_wide)
  d$z_gene_density <- zscore(d$gene_density)
  rownames(d) <- d$tip

  fit <- phylolm::phylolm(
    z_response ~ kingdom + z_log10_genome_size + z_TE_total_pct + z_pi_genome_wide+z_gene_density,
    data = d,
    phy = po$tree,
    model = "lambda",
    REML = FALSE
  )

  beta <- fit$coefficients
  se <- sqrt(diag(fit$vcov))
  z <- beta / se
  p <- 2 * pnorm(abs(z), lower.tail = FALSE)

  tibble(
    root_age_Ma = root_age_Ma,
    response = response,
    n = fit$n,
    lambda = as.numeric(fit$optpar),
    term = names(beta),
    estimate = as.numeric(beta),
    std_error = as.numeric(se),
    ci_low = as.numeric(beta - 1.96 * se),
    ci_high = as.numeric(beta + 1.96 * se),
    P = as.numeric(p),
    AIC = as.numeric(fit$aic),
    R2 = as.numeric(fit$r.squared)
  )
}

# ------------------- phylogenetically stratified resampling -------------------
make_phylo_strata <- function(tree, tips, n_strata) {
  tr <- ape::keep.tip(tree, intersect(tree$tip.label, tips))
  if (ape::Ntip(tr) < n_strata) {
    stop("Requested n_strata exceeds available matched taxa.", call. = FALSE)
  }
  D <- ape::cophenetic.phylo(tr)
  hc <- hclust(as.dist(D), method = "average")
  cl <- cutree(hc, k = n_strata)
  tibble(tip = names(cl), stratum = as.integer(cl))
}

run_phylo_resampling <- function(animal_tree, plant_tree,
                                 animal_dat, plant_dat,
                                 responses,
                                 strata_values = c(50, 75, 100, 125, 150),
                                 B = 1000L) {
  res <- vector("list", length(strata_values) * length(responses))
  idx <- 1L

  for (S in strata_values) {
    a_str <- make_phylo_strata(animal_tree, animal_dat$tip, S)
    p_str <- make_phylo_strata(plant_tree, plant_dat$tip, S)

    a <- animal_dat %>% inner_join(a_str, by = "tip")
    p <- plant_dat %>% inner_join(p_str, by = "tip")

    for (resp in responses) {
      keep_a <- a %>% filter(is.finite(.data[[resp]]))
      keep_p <- p %>% filter(is.finite(.data[[resp]]))

      a_groups <- split(keep_a, keep_a$stratum)
      p_groups <- split(keep_p, keep_p$stratum)

      # Keep only non-empty strata.
      a_groups <- a_groups[vapply(a_groups, nrow, integer(1)) > 0]
      p_groups <- p_groups[vapply(p_groups, nrow, integer(1)) > 0]

      if (length(a_groups) < 10 || length(p_groups) < 10) {
        stop("Too few non-empty phylogenetic strata for response: ", resp, call. = FALSE)
      }

      draws <- replicate(B, {
        av <- vapply(a_groups, function(g) {
          g[[resp]][sample.int(nrow(g), 1L)]
        }, numeric(1))
        pv <- vapply(p_groups, function(g) {
          g[[resp]][sample.int(nrow(g), 1L)]
        }, numeric(1))
        mean(av, na.rm = TRUE) - mean(pv, na.rm = TRUE)
      })

      res[[idx]] <- tibble(
        n_strata_requested = S,
        response = resp,
        B = B,
        mean_difference_Metazoa_minus_Plant = mean(draws),
        median_difference = median(draws),
        ci_low = unname(quantile(draws, 0.025)),
        ci_high = unname(quantile(draws, 0.975)),
        prop_gt_zero = mean(draws > 0),
        prop_lt_zero = mean(draws < 0)
      )
      idx <- idx + 1L
    }
  }
  bind_rows(res)
}

# ------------------------ mapped-vs-unmapped audit -----------------------------
mapped_bias_test <- function(pg, mapped_rows, kingdom_name, variables) {
  d <- pg %>%
    filter(.data$kingdom == kingdom_name) %>%
    mutate(mapping_status = if_else(.data$.row_id %in% mapped_rows,
                                    "Mapped", "Unmapped"))

  bind_rows(lapply(variables, function(v) {
    x <- d %>% filter(is.finite(.data[[v]]))
    n_m <- sum(x$mapping_status == "Mapped")
    n_u <- sum(x$mapping_status == "Unmapped")
    if (n_m < 2 || n_u < 2) {
      return(tibble(
        kingdom = kingdom_name, variable = v,
        n_mapped = n_m, n_unmapped = n_u,
        mapped_median = ifelse(n_m > 0, median(x[[v]][x$mapping_status == "Mapped"]), NA),
        unmapped_median = ifelse(n_u > 0, median(x[[v]][x$mapping_status == "Unmapped"]), NA),
        Wilcoxon_P = NA_real_
      ))
    }
    wt <- wilcox.test(x[[v]] ~ x$mapping_status, exact = FALSE)
    tibble(
      kingdom = kingdom_name, variable = v,
      n_mapped = n_m, n_unmapped = n_u,
      mapped_median = median(x[[v]][x$mapping_status == "Mapped"]),
      unmapped_median = median(x[[v]][x$mapping_status == "Unmapped"]),
      Wilcoxon_P = wt$p.value
    )
  }))
}


# =============================================================================
# READ DATA
# =============================================================================
pg <- readr::read_csv(FILES$per_genome, show_col_types = FALSE, progress = FALSE)
pi_tbl <- readr::read_csv(FILES$per_inversion, show_col_types = FALSE, progress = FALSE)
# Per-genome inversion-divergence summaries.
inv_summary <- pi_tbl %>%
  group_by(sample_id) %>%
  summarise(
    n_inv_pi = sum(is.finite(pi_ratio_inv_genome)),
    median_pi_ratio_inv_genome = safe_median(pi_ratio_inv_genome),
    frac_pi_ratio_inv_genome_gt1 = safe_fraction_gt(pi_ratio_inv_genome, 1),
    median_pi_ratio_inv_flank = safe_median(pi_ratio_inv_flank),
    .groups = "drop"
  )
pg <- pg %>%
  mutate(
    .row_id = row_number(),
    species_tip = gsub("[[:space:]]+", "_", trimws(Species)),
    prop_large = if_else(is.finite(inv_count_total) & inv_count_total > 0,
                         inv_count_large_gt1Mb / inv_count_total, NA_real_),
    log10_inv_count_per_Gb = if_else(
      is.finite(inv_count_per_Gb) & inv_count_per_Gb > 0,
      log10(inv_count_per_Gb), NA_real_
    ),
    log10_genome_size = if_else(
      is.finite(asm_total_len_bp) & asm_total_len_bp > 0,
      log10(asm_total_len_bp), NA_real_
    ),
    logit_prop_large = if_else(
      is.finite(prop_large),
      qlogis(clamp01(prop_large)), NA_real_
    )
  ) %>%
  left_join(inv_summary, by = "sample_id") %>%
  mutate(
    log10_median_pi_ratio = if_else(
      is.finite(median_pi_ratio_inv_genome) & median_pi_ratio_inv_genome > 0,
      log10(median_pi_ratio_inv_genome), NA_real_
    )
  )
# =============================================================================
# READ TREES
# =============================================================================
animal_tree <- read_tree_strict(FILES$animal, "animal")
plant_tree <- read_tree_strict(FILES$plant, "plant")

# =============================================================================
# TREE ↔ GENOME MAPPING
# =============================================================================
animal_map <- map_tree_to_genomes(animal_tree, pg, "Metazoa")
plant_map <- map_tree_to_genomes(plant_tree, pg, "Viridiplantae")

animal_dat <- animal_map$mapped
plant_dat <- plant_map$mapped

coverage_tbl <- tibble(
  kingdom = c("Metazoa", "Viridiplantae"),
  empirical_genomes = c(sum(pg$kingdom == "Metazoa"),
                        sum(pg$kingdom == "Viridiplantae")),
  tree_tips = c(ape::Ntip(animal_tree), ape::Ntip(plant_tree)),
  mapped_genomes = c(nrow(animal_dat), nrow(plant_dat)),
  mapping_coverage_empirical = mapped_genomes / empirical_genomes,
  unmapped_empirical_genomes = empirical_genomes - mapped_genomes,
  unmapped_tree_tips = c(length(animal_map$unmapped_tree_tips),
                         length(plant_map$unmapped_tree_tips))
)
readr::write_csv(coverage_tbl, file.path(OUTDIR, "02_tree_mapping_coverage.csv"))

mapping_detail <- bind_rows(
  animal_dat %>% transmute(kingdom = "Metazoa", tip, sample_id, Species, order,
                           major_clade, mapping_method),
  plant_dat %>% transmute(kingdom = "Viridiplantae", tip, sample_id, Species, order,
                          major_clade, mapping_method)
)
readr::write_csv(mapping_detail, file.path(OUTDIR, "03_tree_mapping_detail.csv"))

mapping_bias <- bind_rows(
  mapped_bias_test(
    pg, animal_dat$row_id, "Metazoa",
    c("log10_inv_count_per_Gb", "median_log10_inv_length",
      "log10_genome_size", "TE_total_pct", "gene_density", "pi_genome_wide")
  ),
  mapped_bias_test(
    pg, plant_dat$row_id, "Viridiplantae",
    c("log10_inv_count_per_Gb", "median_log10_inv_length",
      "log10_genome_size", "TE_total_pct", "gene_density", "pi_genome_wide")
  )
) %>%
  group_by(kingdom) %>%
  mutate(Wilcoxon_FDR = p.adjust(Wilcoxon_P, method = "BH")) %>%
  ungroup()
readr::write_csv(mapping_bias, file.path(OUTDIR, "04_mapping_bias_tests.csv"))


# =============================================================================
# 1. PHYLOGENETIC SIGNAL
# =============================================================================
CORE_TRAITS <- c(
  "log10_inv_count_per_Gb",
  "median_log10_inv_length",
  "logit_prop_large",
  "log10_median_pi_ratio",
  "inv_count_total"
)

signal_tbl <- bind_rows(
  lapply(CORE_TRAITS, function(trait)
    run_phylosig_one(animal_tree, animal_dat, trait, "Metazoa")),
  lapply(CORE_TRAITS, function(trait)
    run_phylosig_one(plant_tree, plant_dat, trait, "Viridiplantae"))
) %>%
  group_by(kingdom) %>%
  mutate(
    lambda_FDR = p.adjust(lambda_P, method = "BH"),
    K_FDR = p.adjust(K_P, method = "BH")
  ) %>%
  ungroup()

readr::write_csv(signal_tbl, file.path(OUTDIR, "05_phylogenetic_signal_FINAL.csv"))

# =============================================================================
# 2. PHYLOGENETIC VARIANCE FRACTION
# =============================================================================
h2_tbl <- bind_rows(
  lapply(CORE_TRAITS, function(trait) {
    profile_h2(animal_tree, animal_dat, trait) %>% mutate(kingdom = "Metazoa")
  }),
  lapply(CORE_TRAITS, function(trait) {
    profile_h2(plant_tree, plant_dat, trait) %>% mutate(kingdom = "Viridiplantae")
  })
) %>%
  select(kingdom, everything())

readr::write_csv(h2_tbl, file.path(OUTDIR, "06_phylogenetic_variance_fraction.csv"))

# =============================================================================
# 3. WITHIN-KINGDOM PGLS
# =============================================================================
PGLS_RESPONSES <- c("log10_inv_count_per_Gb", "median_log10_inv_length")

pgls_primary <- list()
pgls_sens <- list()
k <- 1L

for (king in c("Metazoa", "Viridiplantae")) {
  tr <- if (king == "Metazoa") animal_tree else plant_tree
  dat <- if (king == "Metazoa") animal_dat else plant_dat

  for (resp in PGLS_RESPONSES) {
    pgls_primary[[k]] <- run_pgls(
      tr, dat, resp, king,
      include_gene_density = FALSE,
      model_label = "Primary: genome size + TE + pi"
    )
    pgls_sens[[k]] <- run_pgls(
      tr, dat, resp, king,
      include_gene_density = TRUE,
      model_label = "Sensitivity: + gene density"
    )
    k <- k + 1L
  }
}

primary_coef <- bind_rows(lapply(pgls_primary, `[[`, "coefficients"))
primary_models <- bind_rows(lapply(pgls_primary, `[[`, "model"))
sens_coef <- bind_rows(lapply(pgls_sens, `[[`, "coefficients"))
sens_models <- bind_rows(lapply(pgls_sens, `[[`, "model"))

readr::write_csv(primary_coef, file.path(OUTDIR, "07_PGLS_primary_coefficients.csv"))
readr::write_csv(primary_models, file.path(OUTDIR, "08_PGLS_primary_models.csv"))
readr::write_csv(sens_coef, file.path(OUTDIR, "09_PGLS_gene_density_sensitivity_coefficients.csv"))
readr::write_csv(sens_models, file.path(OUTDIR, "10_PGLS_gene_density_sensitivity_models.csv"))

# Predictor correlations + VIF.
predictor_corr <- bind_rows(lapply(c("Metazoa", "Viridiplantae"), function(king) {
  d <- if (king == "Metazoa") animal_dat else plant_dat
  vars <- c("log10_genome_size", "TE_total_pct", "gene_density", "pi_genome_wide")
  cm <- cor(d[, vars], use = "pairwise.complete.obs", method = "pearson")
  as.data.frame(as.table(cm)) %>%
    as_tibble() %>%
    rename(var1 = Var1, var2 = Var2, r = Freq) %>%
    mutate(kingdom = king, .before = 1)
}))
readr::write_csv(predictor_corr, file.path(OUTDIR, "11_predictor_correlations.csv"))

vif_tbl <- bind_rows(lapply(c("Metazoa", "Viridiplantae"), function(king) {
  d <- if (king == "Metazoa") animal_dat else plant_dat
  bind_rows(lapply(PGLS_RESPONSES, function(resp) {
    bind_rows(
      calc_vif(d,
               c("log10_genome_size", "TE_total_pct", "pi_genome_wide"),
               king, resp, "Primary"),
      calc_vif(d,
               c("log10_genome_size", "TE_total_pct", "pi_genome_wide", "gene_density"),
               king, resp, "Sensitivity + gene density")
    )
  }))
}))
readr::write_csv(vif_tbl, file.path(OUTDIR, "12_VIF_diagnostics.csv"))

# =============================================================================
# 4. MAIN-TEXT CORRELATIONS REPLACEMENT: PHYLOGENETICALLY CORRECTED
# =============================================================================
# Fig.2c: genome-wide pi vs inversion count/Gb
# Fig.2f: genome-wide pi vs median inversion length
main_corr <- bind_rows(
  run_simple_pgls(animal_tree, animal_dat,
                  "log10_inv_count_per_Gb", "pi_genome_wide", "Metazoa"),
  run_simple_pgls(plant_tree, plant_dat,
                  "log10_inv_count_per_Gb", "pi_genome_wide", "Viridiplantae"),
  run_simple_pgls(animal_tree, animal_dat,
                  "median_log10_inv_length", "pi_genome_wide", "Metazoa"),
  run_simple_pgls(plant_tree, plant_dat,
                  "median_log10_inv_length", "pi_genome_wide", "Viridiplantae" ),
  run_simple_pgls(animal_tree, animal_dat,
                  "inv_count_total", "pi_genome_wide", "Metazoa"),
  run_simple_pgls(plant_tree, plant_dat,
                  "inv_count_total", "pi_genome_wide", "Viridiplantae")
  )
readr::write_csv(main_corr, file.path(OUTDIR, "13_main_Fig2_phylogenetic_correlations.csv"))

# =============================================================================
# 5. PHYLOGENETICALLY STRATIFIED RESAMPLING
# =============================================================================
resampling_tbl <- run_phylo_resampling(
  animal_tree, plant_tree,
  animal_dat, plant_dat,
  responses = PGLS_RESPONSES,
  strata_values = c(50, 75, 100, 125, 150),
  B = 1000L
)
readr::write_csv(resampling_tbl,
                 file.path(OUTDIR, "14_phylogenetic_resampling_multiresolution.csv"))


# =============================================================================
# FIGURE EXTENDED DATA — PHYLOGENETIC ROBUSTNESS
# =============================================================================
sig_long <- signal_tbl %>%
  select(kingdom, trait, lambda, K) %>%
  pivot_longer(c(lambda, K), names_to = "metric", values_to = "estimate") %>%
  mutate(
    trait_label = factor(TRAIT_LABELS[trait], levels = rev(unname(TRAIT_LABELS))),
    metric = recode(metric, lambda = "Pagel's lambda", K = "Blomberg's K")
  )
  

pA <- ggplot(sig_long,
             aes(x = estimate, y = trait_label, color = kingdom, shape = kingdom)) +
  geom_vline(xintercept = 0, linetype = 3, linewidth = 0.35, color = "grey60") +
  geom_point(size = 2.1, stroke = 0.4) +
  facet_wrap(~metric, nrow = 1, scales = "free_x") +
  scale_color_manual(values = KINGDOM_COLS) +
  labs(x = "Phylogenetic signal", y = NULL) +
  theme_nature() +
  theme( legend.position = "none")

h2_plot <- h2_tbl %>%
  mutate(trait_label = factor(TRAIT_LABELS[trait],
                              levels = rev(unname(TRAIT_LABELS))))

pB <- ggplot(h2_plot,
             aes(x = H2_phy, y = trait_label, color = kingdom, shape = kingdom)) +
  geom_vline(xintercept = 0, linetype = 3, linewidth = 0.35, color = "grey60") +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high),                               
                 height = 0.15, linewidth = 0.45) +
  geom_point(size = 2.1, stroke = 0.4) +
  scale_color_manual(values = KINGDOM_COLS) +
  coord_cartesian(xlim = c(0, 1)) +
  labs(x = expression(H[phy]^2 ~ "(phylogenetic variance fraction)"), y = NULL) +
  theme_nature() +
  theme(
    axis.text.y = element_blank(), 
    axis.ticks.y = element_blank(),
    legend.position = c(0.90, 0.15), 
    legend.title = element_blank(),
    legend.text = element_text(size = 6) ,
    legend.key.size = grid::unit(3, "mm")
  )


coef_primary_plot <- sens_coef %>%
  filter(term != "(Intercept)") %>%
  mutate(
    term_label = factor(PRED_LABELS[term],
                        levels = rev(c("Genome size", "TE content", "Genome-wide pi","Gene density"))),
    response_label = recode(
      response,
      log10_inv_count_per_Gb = "Inversion count per Gb",
      median_log10_inv_length = "Median inversion length"
    )
  )

pC <- ggplot(coef_primary_plot,
             aes(x = estimate, y = term_label, color = kingdom, shape = kingdom)) +
  geom_vline(xintercept = 0, linetype = 2, linewidth = 0.35, color = "grey50") +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high),
                 height = 0.15, linewidth = 0.45,
                 position = position_dodge(width = 0.35)) +
  geom_point(size = 2.0, stroke = 0.4,
             position = position_dodge(width = 0.35)) +
  facet_wrap(~response_label, nrow = 1, scales = "free_x") +
  scale_color_manual(values = KINGDOM_COLS) +
  labs(x = "Standardized PGLS coefficient", y = NULL) +
  theme_nature() +
  theme(
    legend.position = c(0.95, 0.15), 
    legend.title = element_blank(),
    legend.text = element_text(size = 6) ,
    legend.key.size = grid::unit(3, "mm")
  )


resamp_plot <- resampling_tbl %>%
  mutate(
    response_label = recode(
      response,
      log10_inv_count_per_Gb = "Inversion count per Gb",
      median_log10_inv_length = "Median inversion length"
    ),
    n_strata_requested = factor(n_strata_requested,
                                levels = sort(unique(n_strata_requested)))
  )

pD <- ggplot(resamp_plot,
             aes(x = mean_difference_Metazoa_minus_Plant,
                 y = n_strata_requested)) +
  geom_vline(xintercept = 0, linetype = 2, linewidth = 0.35, color = "grey45") +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high),
                 height = 0.16, linewidth = 0.5, color = "black") +
  geom_point(size = 2.0, color = "black") +
  facet_wrap(~response_label, nrow = 1, scales = "free_x") +
  labs(x = "Metazoa - Viridiplantae", y = "Phylogenetic strata") +
  theme_nature()

root_plot <- cross_root %>%
  filter(term == "kingdomViridiplantae") %>%
  mutate(
    response_label = recode(
      response,
      log10_inv_count_per_Gb = "Inversion count per Gb",
      median_log10_inv_length = "Median inversion length"
    ),
    root_age_Ma = factor(root_age_Ma, levels = c(1200, 1600, 2000),
                         labels = c("1200 Ma", "1600 Ma", "2000 Ma"))
  )



ext_fig  <- ((pA | pB) / (pC | pD)) +
  plot_annotation(tag_levels = "a",
    theme = theme(plot.tag = element_text(face = "bold")))


ggsave(file.path(OUTDIR, "Extended_Data_Fig_phylogenetic_robustness_FINAL.pdf"),
       ext_fig, width = 180, height = 155, units = "mm", device = "pdf")
ggsave(file.path(OUTDIR, "Extended_Data_Fig_phylogenetic_robustness_FINAL.png"),
       ext_fig, width = 180, height = 155, units = "mm", dpi = 300)

