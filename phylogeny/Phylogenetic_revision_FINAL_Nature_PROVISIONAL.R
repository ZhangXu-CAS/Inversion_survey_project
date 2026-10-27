# PROVISIONAL ADAPTATION FOR 2026-10-23 AUTHOR-CONFIRMED DATA (1,174 GENOMES).
# Historical script imported from prior project files; the only intended modifications
# here are the two hardcoded cohort-size preflight assertions (1180=>1174, 187=>181).
# NOT yet executed on final datasets or final dated trees; do not cite numeric results
# without re-running and validating the tree-tip mapping and output manifests.

#!/usr/bin/env Rscript

# =============================================================================
# PHYLOGENETIC ROBUSTNESS — FINAL NATURE REVISION
# Project: Contrasting regimes of chromosomal inversion evolution across plants
#          and animals
#
# Purpose:
#   1) Quantify phylogenetic signal in core inversion traits
#   2) Estimate a phylogenetic variance fraction (H2_phy) by ML variance partition
#   3) Run within-kingdom PGLS controlling genome size, TE content and pi_genome
#   4) Add gene density only as a sensitivity covariate, with VIF diagnostics
#   5) Repeat phylogenetically stratified equal-N resampling across 5 resolutions
#   6) Run cross-kingdom PGLS across 1200/1600/2000 Ma deep-root assumptions
#   7) Produce PGLS replacements for the two main-text OLS correlations
#   8) Quantify phylogeny-mapping attrition and mapped-vs-unmapped bias
#   9) Generate publication-grade Extended Data and supplementary figures
#  10) Generate circular phylogenetic trees with major clade/order/genomic rings
#
# This script contains >15 hard preflight checks and fails early on inconsistent
# inputs. It does not silently repair mismatches.
#
# Run:
#   Rscript Phylogenetic_revision_FINAL_Nature.R
#
# Required packages:
#   ape, phytools, phylolm, dplyr, readr, tidyr, ggplot2, patchwork, scales
# =============================================================================

options(stringsAsFactors = FALSE, warn = 1)
set.seed(20260922)

# ------------------------------ paths -----------------------------------------
get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  hit <- grep("^--file=", args, value = TRUE)
  if (length(hit) == 0) return(normalizePath(getwd(), mustWork = TRUE))
  normalizePath(dirname(sub("^--file=", "", hit[[1]])), mustWork = TRUE)
}
ROOT <- get_script_dir()
setwd(ROOT)

pick_one <- function(candidates, label) {
  hit <- candidates[file.exists(candidates)]
  if (length(hit) == 0) {
    stop("Missing required input [", label, "]. Tried:\n  ",
         paste(candidates, collapse = "\n  "), call. = FALSE)
  }
  normalizePath(hit[[1]], mustWork = TRUE)
}

FILES <- list(
  per_genome = pick_one(c(
    file.path(ROOT, "per_genome.csv"),
    file.path(ROOT, "per_genome(2).csv")
  ), "per_genome"),

  per_inversion = pick_one(c(
    file.path(ROOT, "per_inversion.csv"),
    file.path(ROOT, "per_inversion(3).csv")
  ), "per_inversion"),

  animal = pick_one(c(
    file.path(ROOT, "animal_mcmctree_augmented_verified.nwk"),
    file.path(ROOT, "animal_mcmctree_augmented_verified(1).nwk")
  ), "augmented animal tree"),

  plant = pick_one(c(
    file.path(ROOT, "plant_mcmctree_sanitized.nwk")
  ), "sanitized plant tree"),

  combined_1200 = pick_one(c(
    file.path(ROOT, "combined_plant_animal_root_1200Ma.nwk")
  ), "combined 1200 Ma tree"),

  combined_1600 = pick_one(c(
    file.path(ROOT, "combined_plant_animal_root_1600Ma.nwk")
  ), "combined 1600 Ma tree"),

  combined_2000 = pick_one(c(
    file.path(ROOT, "combined_plant_animal_root_2000Ma.nwk")
  ), "combined 2000 Ma tree")
)

OUTDIR <- file.path(ROOT, "phylogenetic_revision_FINAL")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

# ------------------------------ packages --------------------------------------
REQ_PKGS <- c("ape", "phytools", "phylolm", "dplyr", "readr", "tidyr",
              "ggplot2", "patchwork", "scales")
missing_pkgs <- REQ_PKGS[!vapply(REQ_PKGS, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop(
    "Missing R package(s): ", paste(missing_pkgs, collapse = ", "),
    "\nInstall them before running. Example:\n",
    "install.packages(c(", paste(sprintf('"%s"', missing_pkgs), collapse = ", "), "))",
    call. = FALSE
  )
}
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
})

# ------------------------------ style -----------------------------------------
KINGDOM_COLS <- c(Metazoa = "#E1812C", Viridiplantae = "#00A087")
TRAIT_LABELS <- c(
  log10_inv_count_per_Gb = "Inversion count per Gb",
  median_log10_inv_length = "Median inversion length",
  logit_prop_large = "Proportion >1 Mb",
  log10_median_pi_ratio = "Median pi_inv / pi_genome"
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
                "log10_genome_size", "TE_total_pct", "pi_genome_wide")

  all_dat <- bind_rows(
    animal_dat %>% select(all_of(use_cols)),
    plant_dat %>% select(all_of(use_cols))
  ) %>%
    filter(
      is.finite(.data[[response]]),
      is.finite(.data$log10_genome_size),
      is.finite(.data$TE_total_pct),
      is.finite(.data$pi_genome_wide)
    )

  po <- prune_and_order(tree, all_dat)
  d <- po$data
  d$kingdom <- factor(d$kingdom, levels = c("Metazoa", "Viridiplantae"))
  d$z_response <- zscore(d[[response]])
  d$z_log10_genome_size <- zscore(d$log10_genome_size)
  d$z_TE_total_pct <- zscore(d$TE_total_pct)
  d$z_pi_genome_wide <- zscore(d$pi_genome_wide)
  rownames(d) <- d$tip

  fit <- phylolm::phylolm(
    z_response ~ kingdom + z_log10_genome_size + z_TE_total_pct + z_pi_genome_wide,
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

# ----------------------- circular tree visualization ---------------------------
continuous_cols <- function(x, palette = c("#F7F7F7", "#3C5488")) {
  ok <- is.finite(x)
  out <- rep("#D9D9D9", length(x))
  if (sum(ok) <= 1) return(out)
  z <- scales::rescale(x[ok], to = c(0, 1))
  pal <- grDevices::colorRampPalette(palette)(101)
  out[ok] <- pal[pmax(1, pmin(101, 1 + floor(z * 100)))]
  out
}

clade_palette <- function(x) {
  lev <- sort(unique(as.character(x)))
  n <- length(lev)
  cols <- if (n <= 8) {
    grDevices::hcl.colors(n, "Dark 3")
  } else {
    grDevices::hcl.colors(n, "Dynamic")
  }
  setNames(cols, lev)
}

plot_tree_with_rings <- function(tree, dat, kingdom_name, out_pdf) {
  po <- prune_and_order(tree, dat)
  tr <- po$tree
  d <- po$data

  # Circular tree page.
  pdf(out_pdf, width = 10.5, height = 10.5, useDingbats = FALSE)
  oldpar <- par(no.readonly = TRUE)
  on.exit({par(oldpar); dev.off()}, add = TRUE)
  par(mar = c(0.5, 0.5, 2.3, 0.5), xpd = NA)

  ape::plot.phylo(
    tr, type = "fan", show.tip.label = FALSE,
    edge.color = "grey55", edge.width = 0.25,
    no.margin = FALSE, cex = 0.2,
    label.offset = 0
  )

  plot_env <- get(".PlotPhyloEnv", envir = asNamespace("ape"))
  lp <- get("last_plot.phylo", envir = plot_env)
  xx <- lp$xx[seq_len(ape::Ntip(tr))]
  yy <- lp$yy[seq_len(ape::Ntip(tr))]
  rr <- sqrt(xx^2 + yy^2)
  ux <- xx / rr
  uy <- yy / rr
  r0 <- max(rr, na.rm = TRUE)

  # Ring 1: major clade
  clade_cols <- clade_palette(d$major_clade)
  c1 <- clade_cols[as.character(d$major_clade)]

  # Ring 2: order alternating monochrome by order
  ord_levels <- unique(as.character(d$order))
  ord_cols <- setNames(rep(c("#595959", "#BDBDBD"), length.out = length(ord_levels)),
                       ord_levels)
  c2 <- ord_cols[as.character(d$order)]

  # Ring 3: log10 inversion count/Gb
  c3 <- continuous_cols(d$log10_inv_count_per_Gb,
                        palette = c("#F7F7F7", "#7A5195"))

  # Ring 4: median inversion length
  c4 <- continuous_cols(d$median_log10_inv_length,
                        palette = c("#F7F7F7", "#EF5675"))

  # Ring 5: TE content
  c5 <- continuous_cols(d$TE_total_pct,
                        palette = c("#F7F7F7", "#003F5C"))

  ring_r <- r0 * c(1.035, 1.075, 1.115, 1.155, 1.195)
  ring_cols <- list(c1, c2, c3, c4, c5)
  for (j in seq_along(ring_r)) {
    points(ux * ring_r[j], uy * ring_r[j], pch = 15,
           cex = ifelse(ape::Ntip(tr) > 500, 0.31, 0.48),
           col = ring_cols[[j]])
  }

  # Label the most represented orders only, to avoid illegible labels.
  ord_n <- d %>% count(order, sort = TRUE) %>% filter(!is.na(order), order != "")
  top_orders <- ord_n %>% filter(n >= 5) %>% slice_head(n = 14) %>% pull(order)

  for (od in top_orders) {
    ii <- which(d$order == od)
    if (length(ii) == 0) next
    ang <- atan2(mean(uy[ii]), mean(ux[ii]))
    rx <- cos(ang)
    ry <- sin(ang)
    deg <- ang * 180 / pi
    srt <- deg
    adj <- if (cos(ang) >= 0) c(0, 0.5) else c(1, 0.5)
    if (cos(ang) < 0) srt <- srt + 180
    text(rx * r0 * 1.255, ry * r0 * 1.255, labels = od,
         cex = 0.46, srt = srt, adj = adj)
  }

  title(main = paste0(kingdom_name, " phylogeny (n = ", ape::Ntip(tr), " mapped genomes)"),
        line = 0.5, cex.main = 1.0)

  legend("bottomleft",
         legend = c("Major clade", "Order", "Inversion count/Gb",
                    "Median inversion length", "TE content"),
         pch = 15,
         col = c("#444444", "#999999", "#7A5195", "#EF5675", "#003F5C"),
         bty = "n", cex = 0.66, pt.cex = 0.8)

  # Major-clade key
  legend("bottomright",
         legend = names(clade_cols), pch = 15, col = clade_cols,
         bty = "n", cex = ifelse(length(clade_cols) > 10, 0.47, 0.57),
         ncol = ifelse(length(clade_cols) > 9, 2, 1), pt.cex = 0.7)

  invisible(NULL)
}

# =============================================================================
# READ DATA
# =============================================================================
pg <- readr::read_csv(FILES$per_genome, show_col_types = FALSE, progress = FALSE)
pi_tbl <- readr::read_csv(FILES$per_inversion, show_col_types = FALSE, progress = FALSE)

required_pg <- c(
  "sample_id", "Species", "kingdom", "order", "major_clade",
  "inv_count_total", "inv_count_large_gt1Mb", "inv_count_per_Gb",
  "median_log10_inv_length", "TE_total_pct", "gene_density",
  "asm_total_len_bp", "pi_genome_wide",
  "BUSCO_completeness", "BUSCO_duplication", "asm_scaffold_N50_Mb"
)
required_pi <- c("sample_id", "pi_ratio_inv_genome", "pi_ratio_inv_flank")
require_columns(pg, required_pg, "per_genome")
require_columns(pi_tbl, required_pi, "per_inversion")

numeric_pg <- setdiff(required_pg, c("sample_id", "Species", "kingdom", "order", "major_clade"))
for (nm in numeric_pg) pg[[nm]] <- as_num(pg[[nm]])
pi_tbl$pi_ratio_inv_genome <- as_num(pi_tbl$pi_ratio_inv_genome)
pi_tbl$pi_ratio_inv_flank <- as_num(pi_tbl$pi_ratio_inv_flank)

# Per-genome inversion-divergence summaries.
inv_summary <- pi_tbl %>%
  group_by(sample_id) %>%
  summarise(
    n_inv_pi = sum(is.finite(pi_ratio_inv_genome)),    median_pi_ratio_inv_genome = safe_median(pi_ratio_inv_genome),
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

# Free the 130k-row table after aggregation.
rm(pi_tbl)
invisible(gc())

# =============================================================================
# READ TREES
# =============================================================================
animal_tree <- read_tree_strict(FILES$animal, "animal")
plant_tree <- read_tree_strict(FILES$plant, "plant")
combined_trees <- list(
  `1200` = read_tree_strict(FILES$combined_1200, "combined_1200"),
  `1600` = read_tree_strict(FILES$combined_1600, "combined_1600"),
  `2000` = read_tree_strict(FILES$combined_2000, "combined_2000")
)

tree_qc <- bind_rows(
  tree_audit(animal_tree, "Metazoa"),
  tree_audit(plant_tree, "Viridiplantae"),
  tree_audit(combined_trees[["1200"]], "Combined_1200Ma"),
  tree_audit(combined_trees[["1600"]], "Combined_1600Ma"),
  tree_audit(combined_trees[["2000"]], "Combined_2000Ma")
)

# =============================================================================
# PREFLIGHT: 18 HARD CHECKS
# =============================================================================
preflight <- list()
add_check <- function(name, pass, detail = "") {
  preflight[[length(preflight) + 1L]] <<- tibble(
    check = name, pass = isTRUE(pass), detail = as.character(detail)
  )
}

add_check("1. per_genome rows == 1174", nrow(pg) == 1174, nrow(pg))
add_check("2. per_genome sample_id unique", !anyDuplicated(pg$sample_id),
          paste("duplicates =", sum(duplicated(pg$sample_id))))
add_check("3. expected kingdoms present",
          all(c("Metazoa", "Viridiplantae", "Fungi") %in% unique(pg$kingdom)),
          paste(sort(unique(pg$kingdom)), collapse = ", "))
add_check("4. Metazoa rows == 982", sum(pg$kingdom == "Metazoa") == 982,
          sum(pg$kingdom == "Metazoa"))
add_check("5. Viridiplantae rows == 181", sum(pg$kingdom == "Viridiplantae") == 181,
          sum(pg$kingdom == "Viridiplantae"))
add_check("6. Fungi rows == 11", sum(pg$kingdom == "Fungi") == 11,
          sum(pg$kingdom == "Fungi"))
add_check("7. inversion count/Gb positive where finite",
          all(pg$inv_count_per_Gb[is.finite(pg$inv_count_per_Gb)] > 0))
add_check("8. genome size positive where finite",
          all(pg$asm_total_len_bp[is.finite(pg$asm_total_len_bp)] > 0))
add_check("9. prop_large within [0,1]",
          all(pg$prop_large[is.finite(pg$prop_large)] >= 0 &
              pg$prop_large[is.finite(pg$prop_large)] <= 1))
add_check("10. animal tree unique tips", !anyDuplicated(animal_tree$tip.label),
          ape::Ntip(animal_tree))
add_check("11. plant tree unique tips", !anyDuplicated(plant_tree$tip.label),
          ape::Ntip(plant_tree))
add_check("12. animal tree branch lengths nonnegative",
          all(animal_tree$edge.length >= -1e-10))
add_check("13. plant tree branch lengths nonnegative",
          all(plant_tree$edge.length >= -1e-10))
add_check("14. combined root trees have same tip set",
          setequal(combined_trees[["1200"]]$tip.label, combined_trees[["1600"]]$tip.label) &&
          setequal(combined_trees[["1600"]]$tip.label, combined_trees[["2000"]]$tip.label))
add_check("15. 1200 Ma tree has valid branch lengths",
          all(is.finite(combined_trees[["1200"]]$edge.length)) &&
          all(combined_trees[["1200"]]$edge.length >= -1e-10))
add_check("16. 1600 Ma tree has valid branch lengths",
          all(is.finite(combined_trees[["1600"]]$edge.length)) &&
          all(combined_trees[["1600"]]$edge.length >= -1e-10))
add_check("17. 2000 Ma tree has valid branch lengths",
          all(is.finite(combined_trees[["2000"]]$edge.length)) &&
          all(combined_trees[["2000"]]$edge.length >= -1e-10))
add_check("18. no missing kingdom labels", !anyNA(pg$kingdom) && all(pg$kingdom != ""))

preflight_tbl <- bind_rows(preflight)
readr::write_csv(preflight_tbl, file.path(OUTDIR, "00_preflight_checks.csv"))
readr::write_csv(tree_qc, file.path(OUTDIR, "01_tree_QC.csv"))

if (!all(preflight_tbl$pass)) {
  print(preflight_tbl)
  stop("PREFLIGHT FAILED. See 00_preflight_checks.csv. No analysis was run.", call. = FALSE)
}

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

# Additional post-mapping checks.
if (nrow(animal_dat) < 0.85 * sum(pg$kingdom == "Metazoa")) {
  stop("Animal mapping coverage <85%; inspect mapping before analysis.", call. = FALSE)
}
if (nrow(plant_dat) < 0.85 * sum(pg$kingdom == "Viridiplantae")) {
  stop("Plant mapping coverage <85%; inspect mapping before analysis.", call. = FALSE)
}
if (anyDuplicated(animal_dat$tip) || anyDuplicated(plant_dat$tip)) {
  stop("Duplicate mapped tips detected after preflight.", call. = FALSE)
}

# =============================================================================
# 1. PHYLOGENETIC SIGNAL
# =============================================================================
CORE_TRAITS <- c(
  "log10_inv_count_per_Gb",
  "median_log10_inv_length",
  "logit_prop_large",
  "log10_median_pi_ratio"
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
# 3. WITHIN-KINGDOM PGLS: PRIMARY + GENE-DENSITY SENSITIVITY
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
                  "median_log10_inv_length", "pi_genome_wide", "Viridiplantae")
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
# 6. CROSS-KINGDOM ROOT-AGE SENSITIVITY
# =============================================================================
cross_root <- bind_rows(lapply(names(combined_trees), function(age) {
  tr <- combined_trees[[age]]
  bind_rows(lapply(PGLS_RESPONSES, function(resp) {
    run_cross_kingdom_pgls(
      tr, animal_dat, plant_dat,
      response = resp,
      root_age_Ma = as.numeric(age)
    )
  }))
}))
readr::write_csv(cross_root, file.path(OUTDIR, "15_cross_kingdom_root_sensitivity_PGLS.csv"))

# =============================================================================
# 7. ORDER-LEVEL CONTINUUM
# =============================================================================
order_tbl <- pg %>%
  filter(kingdom %in% c("Metazoa", "Viridiplantae")) %>%
  group_by(kingdom, major_clade, order) %>%
  summarise(
    n = n(),
    mean_log10_inv_count_per_Gb = mean(log10_inv_count_per_Gb, na.rm = TRUE),
    mean_median_log10_inv_length = mean(median_log10_inv_length, na.rm = TRUE),
    mean_log10_genome_size = mean(log10_genome_size, na.rm = TRUE),
    mean_TE_total_pct = mean(TE_total_pct, na.rm = TRUE),
    mean_gene_density = mean(gene_density, na.rm = TRUE),
    mean_pi_genome_wide = mean(pi_genome_wide, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(is.finite(mean_log10_inv_count_per_Gb),
         is.finite(mean_median_log10_inv_length))

readr::write_csv(order_tbl, file.path(OUTDIR, "16_order_level_continuum.csv"))

# =============================================================================
# FIGURE A: EXTENDED DATA — PHYLOGENETIC ROBUSTNESS
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
  theme(legend.position = "top") +
  guides(color = guide_legend(nrow = 1), shape = guide_legend(nrow = 1))

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
  theme(legend.position = "none")

coef_primary_plot <- primary_coef %>%
  filter(term != "(Intercept)") %>%
  mutate(
    term_label = factor(PRED_LABELS[term],
                        levels = rev(c("Genome size", "TE content", "Genome-wide pi"))),
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
  theme(legend.position = "top")

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

pE <- ggplot(root_plot,
             aes(x = estimate, y = root_age_Ma)) +
  geom_vline(xintercept = 0, linetype = 2, linewidth = 0.35, color = "grey45") +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high),
                 height = 0.14, linewidth = 0.5, color = "black") +
  geom_point(size = 2.0, shape = 21, fill = "white", stroke = 0.55) +
  facet_wrap(~response_label, nrow = 1, scales = "free_x") +
  labs(x = "Plant coefficient (standardized PGLS)",
       y = "Assumed deep root") +
  theme_nature()

coverage_plot <- coverage_tbl %>%
  mutate(
    pct = 100 * mapping_coverage_empirical,
    label = sprintf("%.1f%%\n%d/%d", pct, mapped_genomes, empirical_genomes)
  )

pF <- ggplot(coverage_plot,
             aes(x = kingdom, y = pct, fill = kingdom)) +
  geom_col(width = 0.58, color = "black", linewidth = 0.35) +
  geom_text(aes(label = label), vjust = -0.28, size = 2.5) +
  scale_fill_manual(values = KINGDOM_COLS) +
  coord_cartesian(ylim = c(0, 105), clip = "off") +
  labs(x = NULL, y = "Tree mapping coverage (%)") +
  theme_nature() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 20, hjust = 1))

ext_fig <- ((pA | pB) / (pC | pF) / (pD | pE)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold"))

ggsave(file.path(OUTDIR, "Extended_Data_Fig_phylogenetic_robustness_FINAL.pdf"),
       ext_fig, width = 183, height = 215, units = "mm", device = "pdf")
ggsave(file.path(OUTDIR, "Extended_Data_Fig_phylogenetic_robustness_FINAL.png"),
       ext_fig, width = 183, height = 215, units = "mm", dpi = 600)

# =============================================================================
# FIGURE B: GENE-DENSITY SENSITIVITY
# =============================================================================
sens_plot_tbl <- bind_rows(primary_coef, sens_coef) %>%
  filter(term != "(Intercept)") %>%
  mutate(
    term_label = PRED_LABELS[term],
    response_label = recode(
      response,
      log10_inv_count_per_Gb = "Inversion count per Gb",
      median_log10_inv_length = "Median inversion length"
    )
  )

pSens <- ggplot(sens_plot_tbl,
                aes(x = estimate, y = term_label, color = kingdom,
                    shape = model, group = interaction(kingdom, model))) +
  geom_vline(xintercept = 0, linetype = 2, linewidth = 0.35, color = "grey50") +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high),
                 height = 0.13, linewidth = 0.4,
                 position = position_dodge(width = 0.5)) +
  geom_point(size = 1.8, stroke = 0.4,
             position = position_dodge(width = 0.5)) +
  facet_grid(kingdom ~ response_label, scales = "free_x") +
  scale_color_manual(values = KINGDOM_COLS) +
  labs(x = "Standardized PGLS coefficient", y = NULL) +
  theme_nature() +
  theme(legend.position = "top")

ggsave(file.path(OUTDIR, "Supplementary_Fig_gene_density_PGLS_sensitivity.pdf"),
       pSens, width = 180, height = 125, units = "mm", device = "pdf")
ggsave(file.path(OUTDIR, "Supplementary_Fig_gene_density_PGLS_sensitivity.png"),
       pSens, width = 180, height = 125, units = "mm", dpi = 600)

# =============================================================================
# FIGURE C: ORDER-LEVEL CONTINUUM
# =============================================================================
label_orders <- order_tbl %>%
  arrange(desc(n)) %>%
  group_by(kingdom) %>%
  slice_head(n = 12) %>%
  ungroup()

pOrder <- ggplot(order_tbl,
                 aes(x = mean_median_log10_inv_length,
                     y = mean_log10_inv_count_per_Gb,
                     size = n, color = kingdom)) +
  geom_point(alpha = 0.72, stroke = 0.25) +
  geom_text(
    data = label_orders,
    aes(label = order),
    size = 2.0, check_overlap = TRUE,
    nudge_y = 0.035, show.legend = FALSE
  ) +
  scale_color_manual(values = KINGDOM_COLS) +
  scale_size_continuous(range = c(1.3, 6.2)) +
  scale_x_continuous(
    name = "Order mean median inversion length",
    breaks = log10(c(1e2, 1e3, 1e4, 1e5, 1e6)),
    labels = c("0.1 kb", "1 kb", "10 kb", "100 kb", "1 Mb")
  ) +
  scale_y_continuous(
    name = "Order mean inversion count per Gb",
    breaks = log10(c(1, 10, 100, 1000)),
    labels = c("1", "10", "100", "1,000")
  ) +
  labs(size = "Genomes/order") +
  theme_nature(base_size = 8) +
  theme(legend.position = "top")

ggsave(file.path(OUTDIR, "Supplementary_Fig_order_level_continuum.pdf"),
       pOrder, width = 150, height = 120, units = "mm", device = "pdf")
ggsave(file.path(OUTDIR, "Supplementary_Fig_order_level_continuum.png"),
       pOrder, width = 150, height = 120, units = "mm", dpi = 600)

# =============================================================================
# FIGURE D: CIRCULAR PHYLOGENIES
# =============================================================================
plot_tree_with_rings(
  animal_tree, animal_dat, "Metazoa",
  file.path(OUTDIR, "Extended_Data_Tree_Metazoa_FINAL.pdf")
)
plot_tree_with_rings(
  plant_tree, plant_dat, "Viridiplantae",
  file.path(OUTDIR, "Extended_Data_Tree_Viridiplantae_FINAL.pdf")
)

# =============================================================================
# COMPACT RESULTS SUMMARY
# =============================================================================
summary_lines <- c(
  "PHYLOGENETIC REVISION — FINAL AUTOMATED SUMMARY",
  paste0("Generated: ", Sys.time()),
  "",
  "Tree mapping coverage:",
  paste(capture.output(print(coverage_tbl)), collapse = "\n"),
  "",
  "Phylogenetic signal:",
  paste(capture.output(print(signal_tbl)), collapse = "\n"),
  "",
  "Phylogenetic variance fraction:",
  paste(capture.output(print(h2_tbl)), collapse = "\n"),
  "",
  "Primary within-kingdom PGLS coefficients:",
  paste(capture.output(print(primary_coef %>% filter(term != "(Intercept)"))), collapse = "\n"),
  "",
  "Gene-density sensitivity PGLS coefficients:",
  paste(capture.output(print(sens_coef %>% filter(term != "(Intercept)"))), collapse = "\n"),
  "",
  "Multi-resolution phylogenetic resampling:",
  paste(capture.output(print(resampling_tbl)), collapse = "\n"),
  "",
  "Cross-kingdom deep-root sensitivity (kingdom coefficient only):",
  paste(capture.output(print(cross_root %>% filter(term == "kingdomViridiplantae"))), collapse = "\n"),
  "",
  "Main-text phylogenetically corrected correlations:",
  paste(capture.output(print(main_corr)), collapse = "\n"),
  "",
  "Maximum VIF by model:",
  paste(capture.output(
    vif_tbl %>% group_by(kingdom, model) %>%
      summarise(max_VIF = max(VIF, na.rm = TRUE), .groups = "drop")
  ), collapse = "\n")
)

writeLines(summary_lines, file.path(OUTDIR, "RESULTS_SUMMARY.txt"))

# Session information for full reproducibility.
writeLines(capture.output(sessionInfo()), file.path(OUTDIR, "SESSION_INFO.txt"))

message("\nSUCCESS.")
message("All analyses completed. Output directory:\n  ", OUTDIR)
message("\nPrimary figure:\n  ",
        file.path(OUTDIR, "Extended_Data_Fig_phylogenetic_robustness_FINAL.pdf"))
message("\nTrees:\n  ",
        file.path(OUTDIR, "Extended_Data_Tree_Metazoa_FINAL.pdf"), "\n  ",
        file.path(OUTDIR, "Extended_Data_Tree_Viridiplantae_FINAL.pdf"))