############################################################
## Figure 4 (Main) + Figure S5 (Supplementary)
############################################################

library(tidyverse)
library(rstatix)
library(patchwork)
library(scales)
library(dplyr)
library(readr)
library(ggplot2)
library(forcats)
## =========================
## 0) User-editable settings
## =========================
IN_PER_INVERSION <- "per_inversion.csv"
IN_PER_GENOME    <- "per_genome.csv"    
IN_PER_GENE   <- "per_gene_kaks_inversion.csv" 

OUT_DIR          <- "."          # change if you want a dedicated output folder
AGE_THRESHOLD    <- 1.0          # pi_ratio_inv_genome <= 1: Young; >1: Old
BOOT_REPS      <- 5000                       # bootstrap replicates for binned 95% CI
BOOT_CONF      <- 0.95                       # confidence level for bootstrap CI



dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

save_pdf <- function(plot, filename, width_mm, height_mm) {
  ggsave(
    filename = file.path(OUT_DIR, filename),
    plot     = plot,
    width    = width_mm,
    height   = height_mm,
    units    = "mm",
    useDingbats = FALSE
  )
}
save_png <- function(plot, filename, width_mm, height_mm, dpi = 300) {
  ggsave(
    filename = file.path(OUT_DIR, filename),
    plot     = plot,
    width    = width_mm,
    height   = height_mm,
    units    = "mm",
    dpi      = dpi
  )
}

## =========================
## 1) Load & harmonize columns
## =========================
per_inv_raw <- readr::read_csv(IN_PER_INVERSION, show_col_types = FALSE)
per_genome_raw <- readr::read_csv(IN_PER_GENOME, show_col_types = FALSE)

add_missing_cols <- function(df, cols) {
  for (cc in cols) if (!cc %in% names(df)) df[[cc]] <- NA
  df
}

## Add expected columns (including NEW ka_inv / ks_inv)
per_inv_raw <- add_missing_cols(
  per_inv_raw,
  c("sample_id","kingdom","inv_len_bp",
    "gene_count_in_inv","repeat_features_count","te_prop_union",
    "gene_density_inv","gene_density_ratio_genome","gene_density_ratio_flank",
    "pi_inv","pi_ratio_inv_flank","pi_ratio_inv_genome",
    "kaks_median","ka_median","ks_median")
)

## Core per-inversion object used for all plots (animals + plants only)
per_inv <- per_inv_raw %>%
  filter(kingdom %in% c("Metazoa", "Viridiplantae")) %>%
  mutate(
    kingdom    = factor(kingdom, levels = c("Metazoa", "Viridiplantae")),
    
    inv_len_bp = as.numeric(inv_len_bp),
    log10_len  = log10(inv_len_bp),
    
    ## age proxy (π_inversion/π_genome)
    pi_ratio = suppressWarnings(as.numeric(pi_ratio_inv_genome)),
    pi_ratio = if_else(!is.finite(pi_ratio) | pi_ratio <= 0, NA_real_, pi_ratio),
    age_class = case_when(
      is.na(pi_ratio)           ~ NA_character_,
      pi_ratio <= AGE_THRESHOLD ~ "Young",
      pi_ratio >  AGE_THRESHOLD ~ "Old"
    ),
    age_class = factor(age_class, levels = c("Young","Old")),

    gene_count_in_inv         = suppressWarnings(as.numeric(gene_count_in_inv)),
    repeat_features_count     = suppressWarnings(as.numeric(repeat_features_count)),
    te_prop_union             = suppressWarnings(as.numeric(te_prop_union)),
    gene_density_inv          = suppressWarnings(as.numeric(gene_density_inv)),
    gene_density_ratio_genome = suppressWarnings(as.numeric(gene_density_ratio_genome)),
    gene_density_ratio_flank  = suppressWarnings(as.numeric(gene_density_ratio_flank)),

    gc_inv          = suppressWarnings(as.numeric(gc_inv)),
    gc_ratio_flank  = suppressWarnings(as.numeric(gc_ratio_flank)),
    gc_ratio_genome = suppressWarnings(as.numeric(gc_ratio_genome)),

    pi_inv               = suppressWarnings(as.numeric(pi_inv)),
    pi_ratio_inv_flank   = suppressWarnings(as.numeric(pi_ratio_inv_flank)),
    pi_ratio_inv_genome  = suppressWarnings(as.numeric(pi_ratio_inv_genome))
  ) %>%
  filter(is.finite(inv_len_bp), inv_len_bp > 0)


## =========================
## 2) Styling
## =========================
kingdom_cols <- c(
  "Metazoa"       = "#E1812C",
  "Viridiplantae" = "#00A087"
)

theme_nature <- function(base_size = 8) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid      = element_blank(),
      axis.text       = element_text(color = "black"),
      axis.title      = element_text(color = "black"),
      legend.title    = element_blank(),
      panel.border    = element_rect(color = "black", linewidth = 0.5),
      axis.ticks      = element_line(color = "black"),
      plot.title      = element_blank(),
      plot.subtitle   = element_blank()
    )
}

## =========================
## 3) Length-binned helper (log10 bp bins)
## =========================
len_breaks <- seq(2, 8, by = 1)  # log10(bp)

inv_scaling <- per_inv %>%
  mutate(
    len_bin_idx = cut(
      log10_len,
      breaks = len_breaks,
      labels = FALSE,
      include.lowest = TRUE,
      right = FALSE
    )
  ) %>%
  filter(!is.na(len_bin_idx)) %>%
  mutate(
    len_bin_left = len_breaks[len_bin_idx],
    log10_genes = if_else(!is.na(gene_count_in_inv) & gene_count_in_inv > 0,
                          log10(gene_count_in_inv), NA_real_)
  )


## =========================
## Bootstrap CI helper (binned means)
## =========================
boot_ci_mean <- function(x, B = BOOT_REPS, conf = BOOT_CONF) {
  x <- x[is.finite(x) & !is.na(x)]
  n <- length(x)
  if (n < 2) return(c(NA_real_, NA_real_))
  alpha <- (1 - conf) / 2
  reps <- replicate(B, mean(sample(x, size = n, replace = TRUE)))
  as.numeric(stats::quantile(reps, probs = c(alpha, 1 - alpha), na.rm = TRUE, names = FALSE))
}

summarise_boot_ci <- function(df, value_col, group_vars) {
  v <- enquo(value_col)
  g <- rlang::syms(group_vars)

  df %>%
    filter(!is.na(!!v), is.finite(!!v)) %>%
    group_by(!!!g) %>%
    summarise(
      n    = n(),
      mean = mean(!!v),
      sd   = sd(!!v),
      se   = sd / sqrt(n),
      boot = list(boot_ci_mean(pull(cur_data(), !!v))),
      .groups = "drop"
    ) %>%
    mutate(
      ci_low  = purrr::map_dbl(boot, 1),
      ci_high = purrr::map_dbl(boot, 2)
    ) %>%
    select(-boot)
}

summarise_by_lenbin <- function(df, value_col) {
  v <- enquo(value_col)
  dat <- df %>% filter(!is.na(!!v), is.finite(!!v))
  
  summary_tab <- summarise_boot_ci(dat, !!v, group_vars = c('kingdom', 'len_bin_left'))
  
  stat_tab <- dat %>%
    group_by(len_bin_left) %>%
    filter(n_distinct(kingdom) == 2) %>%
    wilcox_test(as.formula(paste0(rlang::as_name(v), " ~ kingdom"))) %>%
    adjust_pvalue(method = "BH") %>%
    add_significance("p.adj") %>%
    ungroup()
  
  bin_y <- summary_tab %>%
    group_by(len_bin_left) %>%
    summarise(y = max(ci_high, na.rm = TRUE), .groups = "drop") %>%
    mutate(y = y * 1.06)
  
  stat_tab <- stat_tab %>% left_join(bin_y, by = "len_bin_left")
  list(summary = summary_tab, stats = stat_tab)
}

plot_lenbin_trait <- function(sum_stat_list, ylab, show_legend = FALSE) {
  summary_tab <- sum_stat_list$summary
  stat_tab    <- sum_stat_list$stats
  
  p <- ggplot(summary_tab,
              aes(x = len_bin_left, y = mean, colour = kingdom, group = kingdom)) +
    geom_line(linewidth = 0.7) +
    geom_point(size = 1.4) +
    geom_errorbar(aes(ymin = ci_low, ymax = ci_high),
                  width = 0.08, linewidth = 0.4) +
    geom_text(
      data = stat_tab,
      aes(x = len_bin_left, y = y, label = p.adj.signif),
      inherit.aes = FALSE,
      size = 2.3
    ) +
    scale_color_manual(values = kingdom_cols) +
    scale_x_continuous(
      breaks = len_breaks,  # 2..8
      labels = len_breaks,
      limits = c(min(len_breaks), max(len_breaks)),
      name   = expression("Inversion length (log"[10] ~ "bp)")
    ) +
    labs(y = ylab) +
    theme_nature() +
    theme(legend.position = "none",axis.text.x = element_text(size = 7),
          axis.text.y = element_text(size = 7))
}
## =========================
## 4) Age-binned helper
## =========================
age_breaks <- c(0, 0.3, 0.5, 0.7, 0.9, 1.1, 1.3, Inf)
age_mids   <- c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4)

inv_age_scaling <- per_inv %>%
  filter(!is.na(pi_ratio), pi_ratio > 0, is.finite(pi_ratio)) %>%
  mutate(
    age = pi_ratio,
    age_bin_idx = cut(
      age,
      breaks = age_breaks,
      labels = FALSE,
      include.lowest = TRUE,
      right = FALSE
    )
  ) %>%
  filter(!is.na(age_bin_idx)) %>%
  mutate(age_mid = age_mids[age_bin_idx])

summarise_by_agebin <- function(df, value_col) {
  v <- enquo(value_col)
  dat <- df %>% filter(!is.na(!!v), is.finite(!!v))
  
  summary_tab <- summarise_boot_ci(dat, !!v, group_vars = c('kingdom', 'age_mid'))
  
  stat_tab <- dat %>%
    group_by(age_mid) %>%
    filter(n_distinct(kingdom) == 2) %>%
    wilcox_test(as.formula(paste0(rlang::as_name(v), " ~ kingdom"))) %>%
    adjust_pvalue(method = "BH") %>%
    add_significance("p.adj") %>%
    ungroup()
  
  bin_y <- summary_tab %>%
    group_by(age_mid) %>%
    summarise(y = max(ci_high, na.rm = TRUE), .groups = "drop") %>%
    mutate(y = y * 1.06)
  
  stat_tab <- stat_tab %>% left_join(bin_y, by = "age_mid")
  list(summary = summary_tab, stats = stat_tab)
}

plot_agebin_trait <- function(sum_stat_list, ylab, show_legend = FALSE) {
  summary_tab <- sum_stat_list$summary
  stat_tab    <- sum_stat_list$stats
  
  p <- ggplot(summary_tab,
              aes(x = age_mid, y = mean, colour = kingdom, group = kingdom)) +
    geom_line(linewidth = 0.7) +
    geom_point(size = 1.4) +
    geom_errorbar(aes(ymin = ci_low, ymax = ci_high),
                  width = 0.03, linewidth = 0.4) +
    geom_text(
      data = stat_tab,
      aes(x = age_mid, y = y, label = p.adj.signif),
      inherit.aes = FALSE,
      size = 2.3
    ) +
    scale_color_manual(values = kingdom_cols) +
    scale_x_continuous(
      breaks = age_mids,
      labels = round(age_mids, 1),
      name   = expression(pi[inversion] / pi[genome])
    ) +
    labs(y = ylab) +
    theme_nature() +
    theme(legend.position = "none", axis.text.x = element_text(size = 7),axis.text.y = element_text(size = 7))
  
}
##kingdom × bin  winsorize
winsorize_by_kingdom_bin <- function(df, value_col, kingdom_col = "kingdom",
                                     bin_col, p = 0.999) {
  value_col   <- rlang::ensym(value_col)
  kingdom_col <- rlang::ensym(kingdom_col)
  bin_col     <- rlang::ensym(bin_col)
  
  df %>%
    dplyr::group_by(!!kingdom_col, !!bin_col) %>%
    dplyr::mutate(
      .q_hi = stats::quantile(!!value_col, probs = p, na.rm = TRUE, names = FALSE),
      .q_lo = stats::quantile(!!value_col, probs = 1 - p, na.rm = TRUE, names = FALSE),
      !!value_col := pmin(pmax(!!value_col, .q_lo), .q_hi)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(-.q_hi, -.q_lo)
}

## =========================
## 5) Length-class boxplot helper (reused for multiple panels)
## =========================
make_size_class <- function(df) {
  df %>%
    mutate(
      size_class = case_when(
        inv_len_bp < 1e4                      ~ "Small",
        inv_len_bp >= 1e4 & inv_len_bp <= 1e6 ~ "Medium",
        inv_len_bp > 1e6                      ~ "Large",
        TRUE                                  ~ NA_character_
      ),
      size_class = factor(size_class, levels = c("Small","Medium","Large"))
    ) %>%
    filter(!is.na(size_class))
}

plot_lenclass_box <- function(df, value_col, ylab,
                              hline1 = NULL,
                              ylim_q = c(0.01, 0.99),
                              dodge_w = 0.65,
                              jitter_w = 0.15,
                              alpha_j = 0.12) {
  v <- enquo(value_col)
  dat <- df %>%
    filter(!is.na(inv_len_bp)) %>%
    make_size_class() %>%
    filter(!is.na(!!v), is.finite(!!v))
  
  yq <- quantile(pull(dat, !!v), ylim_q, na.rm = TRUE)
  
  stat_tab <- dat %>%
    group_by(size_class) %>%
    filter(n_distinct(kingdom) == 2) %>%
    wilcox_test(as.formula(paste0(rlang::as_name(v), " ~ kingdom"))) %>%
    adjust_pvalue(method = "BH") %>%
    add_significance("p.adj") %>%
    ungroup() %>%
    mutate(y = as.numeric(yq[2]) * 1.06)
  
  p <- ggplot(dat, aes(x = size_class, y = !!v, fill = kingdom)) +
    geom_jitter(aes(color = kingdom),
                position = position_jitterdodge(jitter.width = jitter_w, dodge.width = dodge_w),
                size = 0.4, alpha = alpha_j, show.legend = FALSE) +
    geom_boxplot(position = position_dodge(width = dodge_w),
                 width = 0.45, linewidth = 0.4, outlier.shape = NA, alpha = 0.85) +
    scale_fill_manual(values = kingdom_cols) +
    scale_color_manual(values = kingdom_cols) +
    coord_cartesian(ylim = c(as.numeric(yq[1]), as.numeric(yq[2]) * 1.12)) +
    labs(x = "Inversion length class", y = ylab) +
    theme_nature() +
    theme(legend.position = "none", axis.text.x = element_text(size = 7)) +
    geom_text(data = stat_tab,
              aes(x = size_class, y = y, label = p.adj.signif),
              inherit.aes = FALSE, size = 2.3)
  
  if (!is.null(hline1)) {
    p <- p + geom_hline(yintercept = hline1, linetype = "dashed", linewidth = 0.35)
  }
  p
}

## =========================
## 6) Build MAIN Figure 4 panels
## =========================

## Fig4a: gene density ratio (inversion/genome) vs length bin
len_gd <- summarise_by_lenbin(inv_scaling %>% filter(!is.na(gene_density_ratio_genome)),
                              gene_density_ratio_genome)
p1 <- plot_lenbin_trait(len_gd, ylab = "Gene density (inversion / genome)") +
  geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.35)

## Fig4b Gene density ratio vs length class
p2 <- plot_lenclass_box(per_inv, gene_density_ratio_genome,
                          ylab = "Gene density (inversion / genome)",
                          hline1 = 1)

## Fig4c Gene density ratio vs age-bin
age_gd <- summarise_by_agebin(inv_age_scaling %>% filter(!is.na(gene_density_ratio_genome)),
                              gene_density_ratio_genome)
p3 <- plot_agebin_trait(age_gd, ylab = "Gene density (inversion / genome)") +
  geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.35)


##Fig4d TE proportion vs length bin
len_te <- summarise_by_lenbin(inv_scaling %>% filter(!is.na(te_prop_union), te_prop_union >= 0),
                              te_prop_union)

p4 <- plot_lenbin_trait(len_te, ylab = "TE proportion")

## Fig4d: pi_inv/pi_flank vs length bin
inv_scaling_w <- inv_scaling %>%
  dplyr::filter(!is.na(pi_ratio_inv_flank), pi_ratio_inv_flank > 0) %>%
  winsorize_by_kingdom_bin(pi_ratio_inv_flank, kingdom_col = "kingdom",
                           bin_col = len_bin_left, p = 0.999)

len_pi <- summarise_by_lenbin(inv_scaling_w,
                              pi_ratio_inv_flank)
p5 <- plot_lenbin_trait(len_pi, ylab = expression(pi[inversion] / pi[flank])) +
  geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.35)

## Fig4F: pi_inv/pi_flank vs TE-bin
te_breaks <- seq(0, 1, by = 0.2)
te_mids   <- (te_breaks[-length(te_breaks)] + te_breaks[-1]) / 2

inv_TEbin <- inv_scaling %>%
  filter(!is.na(te_prop_union), te_prop_union >= 0, te_prop_union <= 1,
         !is.na(pi_ratio_inv_flank), pi_ratio_inv_flank > 0)%>%
mutate(te_prop=te_prop_union,
  te_bin_idx = cut(te_prop, breaks = te_breaks, labels = FALSE,
                   include.lowest = TRUE, right = FALSE)
) %>%
  filter(!is.na(te_bin_idx)) %>%
  mutate(te_mid = te_mids[te_bin_idx])

summarise_by_TEbin <- function(df, value_col) {
  v <- enquo(value_col)
  dat <- df %>% filter(!is.na(!!v), is.finite(!!v))
  
  summary_tab <- summarise_boot_ci(dat, !!v, group_vars = c('kingdom', 'te_mid'))
  
  stat_tab <- dat %>%
    group_by(te_mid) %>%
    filter(n_distinct(kingdom) == 2) %>%
    wilcox_test(as.formula(paste0(rlang::as_name(v), " ~ kingdom"))) %>%
    adjust_pvalue(method = "BH") %>%
    add_significance("p.adj") %>%
    ungroup()
  
  bin_y <- summary_tab %>%
    group_by(te_mid) %>%
    summarise(y = max(ci_high, na.rm = TRUE), .groups = "drop") %>%
    mutate(y = y * 1.06)
  
  stat_tab <- stat_tab %>% left_join(bin_y, by = "te_mid")
  list(summary = summary_tab, stats = stat_tab)
}

plot_TEbin_trait <- function(sum_stat_list, ylab) {
  summary_tab <- sum_stat_list$summary
  stat_tab    <- sum_stat_list$stats
  
  ggplot(summary_tab,
         aes(x = te_mid, y = mean, colour = kingdom, group = kingdom)) +
    geom_line(linewidth = 0.7) +
    geom_point(size = 1.4) +
    geom_errorbar(aes(ymin = ci_low, ymax = ci_high),
                  width = 0.03, linewidth = 0.4) +
    geom_text(
      data = stat_tab,
      aes(x = te_mid, y = y, label = p.adj.signif),
      inherit.aes = FALSE,
      size = 2.3
    ) +
    scale_color_manual(values = kingdom_cols) +
    scale_x_continuous(
      breaks = te_mids,
      labels = percent(te_mids, accuracy = 1),
      name   = "TE proportion in inversion"
    ) +
    labs(y = ylab) +
    theme_nature() +
    theme(legend.position = "none", axis.text.x = element_text(size = 7)) +
    geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.35)
}


inv_TEbin_w<-  inv_TEbin %>%
  winsorize_by_kingdom_bin(pi_ratio_inv_flank, kingdom_col = "kingdom",
                           bin_col = te_mid, p = 0.999) 
te_pi <- summarise_by_TEbin(inv_TEbin_w, pi_ratio_inv_flank)
p6 <- plot_TEbin_trait(te_pi, ylab = expression(pi[inversion] / pi[flank]))

## =========================
## 7) ω outlier panels + gene-number matched null models
## =========================

N_NULL <- 100        # number of null resamples per inversion
SEED_BASE <- 20260326 # reproducibility

size_cols <- c(
  "small"  = "#9ecae1",
  "medium" = "#fdae6b",
  "large"  = "#a1d99b"
)

## -------------------------
## Read data
## -------------------------
gene_omega <- read_csv(IN_PER_GENE, show_col_types = FALSE) %>%
  mutate(
    omega = suppressWarnings(as.numeric(omega)),
    omega = ifelse(is.finite(omega) & omega >= 0, omega, NA_real_)
  ) %>%
  filter(!is.na(sample_id), !is.na(inv_id), !is.na(gene_id), !is.na(omega))

gene_omega <- gene_omega %>%
  mutate(gene_id_old = gene_id,
         gene_id = paste(sample_id, gene_id, sep = "_"))

inv_meta <- read_csv(IN_PER_INVERSION, show_col_types = FALSE) %>%
  mutate(
    inv_len_bp = suppressWarnings(as.numeric(inv_len_bp)),
    kingdom = factor(kingdom, levels = c("Metazoa", "Viridiplantae")),
    size_bin = case_when(
      inv_len_bp < 1e5 ~ "small",
      inv_len_bp < 1e6 ~ "medium",
      TRUE ~ "large"
    ),
    size_bin = factor(size_bin, levels = c("small", "medium", "large"))
  ) %>%
  select(sample_id, inv_id, kingdom, inv_len_bp, size_bin)

gene_join <- gene_omega %>%
  inner_join(inv_meta, by = c("sample_id", "inv_id")) %>%
  filter(!is.na(kingdom), !is.na(size_bin))


## -------------------------
## Observed inversion-level ω summaries
## -------------------------
inv_omega_obs <- gene_join %>%
  group_by(sample_id, inv_id, kingdom, size_bin, inv_len_bp) %>%
  summarise(
    n_genes   = n(),
    omega_p95 = quantile(omega, probs = 0.95, na.rm = TRUE, type = 7),
    omega_p99 = quantile(omega, probs = 0.99, na.rm = TRUE, type = 7),
    omega_max = max(omega, na.rm = TRUE),
    omega_med = median(omega, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    across(c(omega_p95, omega_p99, omega_max, omega_med),
           ~ ifelse(is.finite(.x) & .x >= 0, .x, NA_real_))
  ) %>%
  filter(!is.na(omega_p95), !is.na(omega_p99), n_genes >= 1)
## -------------------------
## Build sample-level gene pools
## Null is drawn from the same sample to preserve sample-specific ω background
## -------------------------
sample_gene_pool <- gene_omega %>%
  group_by(sample_id) %>%
  summarise(
    omega_pool = list(omega[is.finite(omega) & !is.na(omega)]),
    n_pool = sum(is.finite(omega) & !is.na(omega)),
    .groups = "drop"
  )

inv_omega_obs <- inv_omega_obs %>%
  left_join(sample_gene_pool, by = "sample_id") %>%
  filter(n_pool >= n_genes)

## -------------------------
## Helper: null distribution for matched gene number
## -------------------------
calc_null_stats <- function(pool, n_genes, nrep = 100, seed = 1) {
  if (length(pool) < n_genes || n_genes < 1) {
    return(tibble(
      null_mean_p95 = NA_real_,
      null_sd_p95   = NA_real_,
      null_mean_p99 = NA_real_,
      null_sd_p99   = NA_real_
    ))
  }
  
  set.seed(seed)
  
  ## replicate returns matrix if simplify = TRUE and fixed-length numeric output
  p95_vals <- replicate(
    nrep,
    quantile(sample(pool, size = n_genes, replace = FALSE),
             probs = 0.95, na.rm = TRUE, type = 7)
  )
  
  p99_vals <- replicate(
    nrep,
    quantile(sample(pool, size = n_genes, replace = FALSE),
             probs = 0.99, na.rm = TRUE, type = 7)
  )
  
  tibble(
    null_mean_p95 = mean(p95_vals, na.rm = TRUE),
    null_sd_p95   = sd(p95_vals, na.rm = TRUE),
    null_mean_p99 = mean(p99_vals, na.rm = TRUE),
    null_sd_p99   = sd(p99_vals, na.rm = TRUE)
  )
}

## -------------------------
## Run matched null per inversion
## -------------------------

null_tbl <- purrr::pmap_dfr(
  list(inv_omega_obs$omega_pool, inv_omega_obs$n_genes, seq_len(nrow(inv_omega_obs))),
  function(pool, n_genes, idx) {
    calc_null_stats(
      pool = pool,
      n_genes = n_genes,
      nrep = N_NULL,
      seed = SEED_BASE + idx
    )
  }
)

inv_omega <- bind_cols(inv_omega_obs %>% select(-omega_pool, -n_pool), null_tbl) %>%
  mutate(
    ratio_p95 = omega_p95/null_mean_p95,
    ratio_p99 = omega_p99/null_mean_p99,
    z_p95 = ifelse(is.finite(null_sd_p95) & null_sd_p95 > 0,
                   (omega_p95 - null_mean_p95) / null_sd_p95, NA_real_),
    z_p99 = ifelse(is.finite(null_sd_p99) & null_sd_p99 > 0,
                   (omega_p99 - null_mean_p99) / null_sd_p99, NA_real_)
  )

## -------------------------
## Helper functions
## -------------------------
p_to_star <- function(p) {
  dplyr::case_when(
    is.na(p) ~ "",
    p < 1e-3 ~ "***",
    p < 1e-2 ~ "**",
    p < 5e-2 ~ "*",
    TRUE     ~ "ns"
  )
}


safe_density_max <- function(x) {
  x <- x[is.finite(x)]
  if (length(unique(x)) < 2) return(1)
  max(stats::density(x, na.rm = TRUE)$y, na.rm = TRUE)
}

plot_density_kingdom <- function(df, xcol, xlab, xlim = c(0, 3)) {
  dat <- df %>% filter(is.finite(.data[[xcol]]))
  wt <- wilcox.test(dat[[xcol]] ~ dat$kingdom)$p.value
  lab <- paste0(ifelse(wt < 0.001, "p < 0.001***", paste0("p = ", signif(wt, 2), " ", p_to_star(wt))))
  dens_max <- max(stats::density(dat[[xcol]], na.rm = TRUE)$y, na.rm = TRUE)
  
  med_tbl <- dat %>%
    group_by(kingdom) %>%
    summarise(med = median(.data[[xcol]], na.rm = TRUE), .groups = "drop")
  
  ggplot(dat, aes(x = .data[[xcol]], color = kingdom, fill = kingdom)) +
    geom_density(alpha = 0.25, adjust = 1.0, linewidth = 0.7) +
    geom_vline(data = med_tbl, aes(xintercept = med, color = kingdom),
               linewidth = 0.6, linetype = "dashed", show.legend = FALSE) +
    scale_color_manual(values = kingdom_cols) +
    scale_fill_manual(values = kingdom_cols) +
    coord_cartesian(xlim = xlim) +
    annotate("text", x = xlim[2] - 0.02 * diff(xlim), y = dens_max * 0.75,
             label = lab, hjust = 1, vjust = 1, size = 2.6, color = "black") +
    labs(x = xlab, y = "Density") +
    theme_nature() +
    theme(legend.position = c(0.98, 0.98),
          legend.justification = c(1, 1),
          legend.direction = "vertical",  
          legend.text = element_text(size = 6),
          legend.key.size = unit(3, "mm"))
}
plot_density_size <- function(df, xcol, xlab, xlim = c(0, 3)) {
  dat <- df %>% filter(is.finite(.data[[xcol]]))
  
  # KW across size bins (within the filtered data)
  kw <- kruskal.test(dat[[xcol]] ~ dat$size_bin)$p.value
  lab <- paste0(
    ifelse(kw < 0.001, "p < 0.001***", paste0("p = ", signif(kw, 2), " ", p_to_star(kw)))
  )
  
  # medians per size_bin (used for dashed vlines)
  med_size <- dat %>%
    group_by(size_bin) %>%
    summarise(med = median(.data[[xcol]], na.rm = TRUE), .groups = "drop")
  
  # For annotation placement (global; OK when facets share y scale)
  dens_max <- max(stats::density(dat[[xcol]], na.rm = TRUE)$y, na.rm = TRUE)
  
  ggplot(dat, aes(x = .data[[xcol]], color = size_bin, fill = size_bin)) +
    geom_density(alpha = 0.20, adjust = 1.0, linewidth = 0.7) +
    geom_vline(
      data = med_size,
      aes(xintercept = med, color = size_bin),
      linewidth = 0.6, linetype = "dashed", show.legend = FALSE
    ) +
    facet_wrap(~ kingdom) +                  
    coord_cartesian(xlim = xlim) +             
    annotate(
      "text",
      x = xlim[2] - 0.02 * diff(xlim),
      y = dens_max ,
      label = lab, hjust = 1, vjust = 1,
      size = 2.6, color = "black"
    ) +
    labs(x = xlab, y = "Density") +
    theme_nature() +
    theme(
      legend.position = c(0.98, 0.98),
      legend.justification = c(1, 1),
      legend.direction = "vertical",
      legend.text = element_text(size = 6),
      legend.key.size = unit(3, "mm")
    )
}


##p7 & p8 
p7<- plot_density_kingdom(inv_omega, "omega_p99",
                          xlab = expression(paste("Null-adjusted Inversion-level ", omega, " (95th percentile)")),xlim = c(0, 3))

p8 <- plot_density_size(inv_omega, "omega_p99",
                        xlab = expression(paste("Null-adjusted Inversion-level ", omega, " (95th percentile)")),xlim = c(0, 3))

##
## =========================
## 7) Simple Ka/Ks shift panels
##    Replace your current omega panel section with this block
## =========================
IN_KAKS_SHIFT <- "inversion_kaks_shift.csv"
size_cols <- c(
  "small"  = "#9ecae1",
  "medium" = "#fdae6b",
  "large"  = "#a1d99b"
)

p_to_star <- function(p) {
  dplyr::case_when(
    is.na(p) ~ "",
    p < 1e-3 ~ "***",
    p < 1e-2 ~ "**",
    p < 5e-2 ~ "*",
    TRUE     ~ "ns"
  )
}

kaks_shift <- readr::read_tsv(IN_KAKS_SHIFT, show_col_types = FALSE) %>%
  mutate(
    inv_kaks_median        = suppressWarnings(as.numeric(inv_kaks_median)),
    inv_kaks_p90           = suppressWarnings(as.numeric(inv_kaks_p90)),
    bg_kaks_median         = suppressWarnings(as.numeric(bg_kaks_median)),
    bg_kaks_p90           = suppressWarnings(as.numeric(bg_kaks_p90)),
    kaks_shift_ratio_median = suppressWarnings(as.numeric(kaks_shift_ratio_median)),
    kaks_shift_ratio_p90    = suppressWarnings(as.numeric(kaks_shift_ratio_p90)),
    kaks_shift_delta_median = suppressWarnings(as.numeric(kaks_shift_delta_median)),
    kaks_shift_delta_p90    = suppressWarnings(as.numeric(kaks_shift_delta_p90)),
    n_genes_found           = suppressWarnings(as.numeric(n_genes_found))
  )

inv_meta_shift <- per_inv_raw %>%
  mutate(
    inv_len_bp = suppressWarnings(as.numeric(inv_len_bp)),
    kingdom = factor(kingdom, levels = c("Metazoa", "Viridiplantae")),
    size_bin = case_when(
      inv_len_bp < 1e5 ~ "small",
      inv_len_bp < 1e6 ~ "medium",
      TRUE ~ "large"
    ),
    size_bin = factor(size_bin, levels = c("small", "medium", "large"))
  ) %>%
  select(sample_id, inv_id, kingdom, inv_len_bp, size_bin)

shift_df <- kaks_shift %>%
  inner_join(inv_meta_shift, by = c("sample_id", "inv_id")) %>%
  filter(
    kingdom %in% c("Metazoa", "Viridiplantae"),
    status == "ok",
    !is.na(size_bin),
    !is.na(n_genes_found), n_genes_found >= 1
  )

plot_density_kingdom_ratio <- function(df, xcol, xlab, xlim = c(0, 3), vline = 1) {
  dat <- df %>% filter(is.finite(.data[[xcol]]), !is.na(.data[[xcol]]))
  wt <- wilcox.test(dat[[xcol]] ~ dat$kingdom)$p.value
  lab <- ifelse(wt < 0.001, "p < 0.001***", paste0("p = ", signif(wt, 2), " ", p_to_star(wt)))
  dens_max <- max(stats::density(dat[[xcol]], na.rm = TRUE)$y, na.rm = TRUE)
  
  med_tbl <- dat %>%
    group_by(kingdom) %>%
    summarise(med = median(.data[[xcol]], na.rm = TRUE), .groups = "drop")
  
  ggplot(dat, aes(x = .data[[xcol]], color = kingdom, fill = kingdom)) +
    geom_density(alpha = 0.25, adjust = 1.0, linewidth = 0.7) +
    geom_vline(data = med_tbl, aes(xintercept = med, color = kingdom),
               linewidth = 0.6, linetype = "dashed", show.legend = FALSE) +
    scale_color_manual(values = kingdom_cols) +
    scale_fill_manual(values = kingdom_cols) +
    coord_cartesian(xlim = xlim) +
    annotate("text", x = xlim[2] - 0.02 * diff(xlim), y = dens_max * 0.75,
             label = lab, hjust = 1, vjust = 1, size = 2.6, color = "black") +
    labs(x = xlab, y = "Density") +
    theme_nature() +
    theme(
      legend.position = c(0.98, 0.98),
      legend.justification = c(1, 1),
      legend.direction = "vertical",
      legend.text = element_text(size = 6),
      legend.key.size = unit(3, "mm")
    )
}
plot_density_size <- function(df, xcol, xlab, xlim = NULL) {
  dat <- df %>% filter(is.finite(.data[[xcol]]), !is.na(size_bin))
  
  kw <- tryCatch(
    kruskal.test(dat[[xcol]] ~ dat$size_bin)$p.value,
    error = function(e) NA_real_
  )
  
  lab <- paste0(
    ifelse(is.na(kw), "NA",
           ifelse(kw < 0.001, "p < 0.001***",
                  paste0("p = ", signif(kw, 2), " ", p_to_star(kw))))
  )
  
  dens_max <- safe_density_max(dat[[xcol]])
  
  med_size <- dat %>%
    group_by(size_bin) %>%
    summarise(med = median(.data[[xcol]], na.rm = TRUE), .groups = "drop")
  
  p <- ggplot(dat, aes(x = .data[[xcol]], color = size_bin, fill = size_bin)) +
    geom_density(alpha = 0.20, adjust = 1.0, linewidth = 0.7) +
    geom_vline(
      data = med_size,
      aes(xintercept = med, color = size_bin),
      linewidth = 0.6, linetype = "dashed", show.legend = FALSE
    ) +
    scale_color_manual(values = size_cols) +
    scale_fill_manual(values = size_cols) +
    facet_wrap(~ kingdom) +
    labs(x = xlab, y = "Density") +
    theme_nature() +
    theme(
      legend.position = c(0.98, 0.98),
      legend.justification = c(1, 1),
      legend.direction = "vertical",
      legend.text = element_text(size = 6),
      legend.key.size = unit(3, "mm")
    )
  
  if (!is.null(xlim)) {
    p <- p + coord_cartesian(xlim = xlim)
    xpos <- xlim[2] - 0.02 * diff(xlim)
  } else {
    xr <- range(dat[[xcol]], na.rm = TRUE)
    xpos <- xr[2] - 0.02 * diff(xr)
  }
  
  p + annotate(
    "text",
    x = xpos, y = dens_max * 2,
    label = lab, hjust = 1, vjust = 1,
    size = 2.6, color = "black"
  )
}

## Main figure panels: use median-based ratio as primary statistic
p7 <- plot_density_kingdom_ratio(
  shift_df,
  "kaks_shift_ratio_p90",
  xlab = expression(paste("Relative ", omega, " of inversion genes (90th percentile)")),
  xlim = c(0, 2)
)

p8<- plot_density_size(
  shift_df,
  "kaks_shift_ratio_p90",
  xlab = expression(paste("Relative ", omega, " of inversion genes (90th percentile)")),
  xlim = c(0, 2)
)

#### Fig4 all
row3 <- p7 + p8 + plot_layout(widths = c(1, 2))
Fig4_all <- (p1 | p2 | p3) /
  (p4 | p5 | p6) /
  row3 +
  plot_annotation(
    tag_levels = "a",
    tag_prefix = "",
    tag_sep    = "",
    theme = theme(plot.tag = element_text(face = "bold", size = 8))
  )
save_pdf(Fig4_all, "Fig4_all.pdf", 180, 180)
save_png(Fig4_all, "Fig4_all.png", 180, 180, dpi = 300)
###Figure S5 panles####

## FigS5A content (genes per inversion; log10) — used in S5 panel a
len_gene <- summarise_by_lenbin(inv_scaling, log10_genes)
pS5_A <- plot_lenbin_trait(
  len_gene,
  ylab = expression("Number of genes per inversion (log"[10] * ")"),
  show_legend = TRUE
)
###pS5_B TE lenclass_box
pS5_B <- plot_lenclass_box(per_inv, te_prop_union, ylab = "TE proportion")

## S5C: pi_inv/pi_flank vs age-bin
inv_age_scaling_w <- inv_age_scaling %>%
  dplyr::filter(!is.na(pi_ratio_inv_flank), pi_ratio_inv_flank > 0) %>%
  winsorize_by_kingdom_bin(pi_ratio_inv_flank, kingdom_col = "kingdom",
                           bin_col = age_mid, p = 0.999)

age_pi <- summarise_by_agebin(inv_age_scaling_w, pi_ratio_inv_flank)

pS5_C <- plot_agebin_trait(age_pi, ylab = expression(pi[inversion] / pi[flank])) +
  geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.35) 

## 

pS5_D <- plot_density_kingdom_ratio(
  shift_df,
  "kaks_shift_ratio_median",
  xlab = expression(paste("Relative ", omega, " of inversion genes")),
  xlim = c(0, 2)
)

pS5_E <- plot_density_size(
  shift_df,
  "kaks_shift_ratio_median",
  xlab = expression(paste("Relative ", omega, " of inversion genes")),
  xlim = c(0, 2)
)


##FigS5_new
row2 <- pS5_D + pS5_E + plot_layout(widths = c(1, 2))

FigS5_new <- (pS5_A |pS5_B | pS5_C)/ row2  +
   plot_annotation(
    tag_levels = "a",
    tag_prefix = "",
    tag_sep    = "",
    theme = theme(plot.tag = element_text(face = "bold", size = 8))
  )

## Taller figure to accommodate 3 rows
save_pdf(FigS5_new, "FigS5_new.pdf", 180, 120)
save_png(FigS5_new, "FigS5_new.png", 180, 120, dpi = 300)

## =========================
## 10) Reproducibility
## =========================
writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))
