
suppressPackageStartupMessages({
  library(tidyverse)
  library(scales)
  library(rstatix)
  library(ggpubr)
  library(patchwork)
})

AGE_THRESHOLD_DEFAULT <- 1.0
BOOT_REPS_DEFAULT <- 5000
BOOT_CONF_DEFAULT <- 0.95
MAX_POINTS_PER_KINGDOM_DEFAULT <- 10000L

load_per_inversion_age <- function(path = IN_PER_INVERSION,
                                   kingdoms = c("Metazoa", "Viridiplantae"),
                                   age_threshold = AGE_THRESHOLD_DEFAULT) {
  per_inv_raw <- readr::read_csv(path, show_col_types = FALSE)
  if (!("pi_ratio_inv_genome" %in% names(per_inv_raw)) && ("pi_ration_inv_genome" %in% names(per_inv_raw))) {
    per_inv_raw <- per_inv_raw %>% rename(pi_ratio_inv_genome = pi_ration_inv_genome)
  }
  add_missing <- c("sample_id","inv_id","gene_count_in_inv","te_prop_union","gene_density_ratio_genome","pi_ratio_inv_flank","pi_ratio_inv_genome")
  for (cc in add_missing) if (!cc %in% names(per_inv_raw)) per_inv_raw[[cc]] <- NA
  per_inv_raw %>%
    filter(kingdom %in% kingdoms, !is.na(inv_len_bp)) %>%
    mutate(
      kingdom = factor(kingdom, levels = kingdoms),
      inv_len_bp = as.numeric(inv_len_bp),
      log10_len = safe_log10(inv_len_bp),
      pi_ratio = suppressWarnings(as.numeric(pi_ratio_inv_genome)),
      pi_ratio = if_else(!is.finite(pi_ratio) | pi_ratio <= 0, NA_real_, pi_ratio),
      age_class = case_when(is.na(pi_ratio) ~ NA_character_, pi_ratio <= age_threshold ~ "Young", pi_ratio > age_threshold ~ "Old"),
      age_class = factor(age_class, levels = c("Young","Old")),
      size_class = case_when(inv_len_bp < 1e4 ~ "Small", inv_len_bp >= 1e4 & inv_len_bp <= 1e6 ~ "Medium", inv_len_bp > 1e6 ~ "Large", TRUE ~ NA_character_),
      size_class = factor(size_class, levels = c("Small","Medium","Large")),
      gene_count_in_inv = suppressWarnings(as.numeric(gene_count_in_inv)),
      log10_genes = if_else(!is.na(gene_count_in_inv) & gene_count_in_inv > 0, safe_log10(gene_count_in_inv), NA_real_),
      te_prop_union = suppressWarnings(as.numeric(te_prop_union)),
      gene_density_ratio_genome = suppressWarnings(as.numeric(gene_density_ratio_genome)),
      pi_ratio_inv_flank = suppressWarnings(as.numeric(pi_ratio_inv_flank))
    ) %>% filter(is.finite(inv_len_bp), inv_len_bp > 0)
}

safe_density_max <- function(x) {
  x <- x[is.finite(x) & !is.na(x)]
  if (length(unique(x)) < 2) return(1)
  max(stats::density(x, na.rm = TRUE)$y, na.rm = TRUE)
}

boot_ci_mean <- function(x, B = BOOT_REPS_DEFAULT, conf = BOOT_CONF_DEFAULT) {
  x <- x[is.finite(x) & !is.na(x)]
  n <- length(x)
  if (n < 2) return(c(NA_real_, NA_real_))
  alpha <- (1 - conf) / 2
  reps <- replicate(B, mean(sample(x, size = n, replace = TRUE)))
  as.numeric(stats::quantile(reps, probs = c(alpha, 1 - alpha), na.rm = TRUE, names = FALSE))
}

summarise_boot_ci <- function(df, value_col, group_vars) {
  v <- enquo(value_col); g <- rlang::syms(group_vars)
  df %>% filter(!is.na(!!v), is.finite(!!v)) %>% group_by(!!!g) %>%
    summarise(n = n(), mean = mean(!!v), sd = sd(!!v), se = sd/sqrt(n), boot = list(boot_ci_mean(pull(cur_data(), !!v))), .groups = "drop") %>%
    mutate(ci_low = purrr::map_dbl(boot, 1), ci_high = purrr::map_dbl(boot, 2)) %>% select(-boot)
}

winsorize_by_group_bin <- function(df, value_col, group_col = "kingdom", bin_col, p = 0.999) {
  value_col <- rlang::ensym(value_col); group_col <- rlang::ensym(group_col); bin_col <- rlang::ensym(bin_col)
  df %>% group_by(!!group_col, !!bin_col) %>%
    mutate(.q_hi = stats::quantile(!!value_col, probs = p, na.rm = TRUE, names = FALSE),
           .q_lo = stats::quantile(!!value_col, probs = 1 - p, na.rm = TRUE, names = FALSE),
           !!value_col := pmin(pmax(!!value_col, .q_lo), .q_hi)) %>%
    ungroup() %>% select(-.q_hi, -.q_lo)
}

summarise_by_lenbin <- function(df, value_col, len_breaks = seq(2, 8, by = 1)) {
  v <- enquo(value_col)
  dat <- df %>% mutate(len_bin_idx = cut(log10_len, breaks = len_breaks, labels = FALSE, include.lowest = TRUE, right = FALSE), len_bin_left = len_breaks[len_bin_idx]) %>%
    filter(!is.na(len_bin_idx), !is.na(!!v), is.finite(!!v))
  summary_tab <- summarise_boot_ci(dat, !!v, group_vars = c('kingdom','len_bin_left'))
  stat_tab <- dat %>% group_by(len_bin_left) %>% filter(n_distinct(kingdom) == 2) %>% wilcox_test(as.formula(paste0(rlang::as_name(v),' ~ kingdom'))) %>% adjust_pvalue(method='BH') %>% add_significance('p.adj') %>% ungroup()
  bin_y <- summary_tab %>% group_by(len_bin_left) %>% summarise(y = max(ci_high, na.rm=TRUE), .groups='drop') %>% mutate(y = y * 1.06)
  list(summary = summary_tab, stats = left_join(stat_tab, bin_y, by='len_bin_left'))
}

plot_lenbin_trait <- function(sum_stat_list, ylab, hline1 = NULL) {
  summary_tab <- sum_stat_list$summary; stat_tab <- sum_stat_list$stats
  p <- ggplot(summary_tab, aes(x = len_bin_left, y = mean, colour = kingdom, group = kingdom)) +
    geom_line(linewidth = 0.7) + geom_point(size = 1.4) + geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.08, linewidth = 0.4) +
    geom_text(data = stat_tab, aes(x = len_bin_left, y = y, label = p.adj.signif), inherit.aes = FALSE, size = 2.3) +
    scale_color_manual(values = kingdom_cols) +
    scale_x_continuous(breaks = seq(2,8,1), labels = seq(2,8,1), limits = c(2,8), name = expression('Inversion length (log'[10] ~ 'bp)')) +
    labs(y = ylab) + theme_nature() + theme(legend.position = 'none', axis.text.x = element_text(size=7), axis.text.y = element_text(size=7))
  if (!is.null(hline1)) p <- p + geom_hline(yintercept = hline1, linetype = 'dashed', linewidth = 0.35)
  p
}

summarise_by_agebin <- function(df, value_col, age_breaks = c(0,0.3,0.5,0.7,0.9,1.1,1.3,Inf), age_mids = c(0.2,0.4,0.6,0.8,1.0,1.2,1.4)) {
  v <- enquo(value_col)
  dat <- df %>% filter(!is.na(pi_ratio), pi_ratio > 0, is.finite(pi_ratio), !is.na(!!v), is.finite(!!v)) %>% mutate(age_mid = age_mids[cut(pi_ratio, breaks = age_breaks, labels = FALSE, include.lowest = TRUE, right = FALSE)]) %>% filter(!is.na(age_mid))
  summary_tab <- summarise_boot_ci(dat, !!v, group_vars = c('kingdom','age_mid'))
  stat_tab <- dat %>% group_by(age_mid) %>% filter(n_distinct(kingdom)==2) %>% wilcox_test(as.formula(paste0(rlang::as_name(v),' ~ kingdom'))) %>% adjust_pvalue(method='BH') %>% add_significance('p.adj') %>% ungroup()
  bin_y <- summary_tab %>% group_by(age_mid) %>% summarise(y = max(ci_high, na.rm=TRUE), .groups='drop') %>% mutate(y = y * 1.06)
  list(summary = summary_tab, stats = left_join(stat_tab, bin_y, by='age_mid'))
}

plot_agebin_trait <- function(sum_stat_list, ylab, hline1 = NULL) {
  summary_tab <- sum_stat_list$summary; stat_tab <- sum_stat_list$stats
  p <- ggplot(summary_tab, aes(x = age_mid, y = mean, colour = kingdom, group = kingdom)) +
    geom_line(linewidth = 0.7) + geom_point(size = 1.4) + geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.03, linewidth = 0.4) +
    geom_text(data = stat_tab, aes(x = age_mid, y = y, label = p.adj.signif), inherit.aes = FALSE, size = 2.3) +
    scale_color_manual(values = kingdom_cols) + scale_x_continuous(breaks = c(0.2,0.4,0.6,0.8,1.0,1.2,1.4), labels = c(0.2,0.4,0.6,0.8,1.0,1.2,1.4), name = expression(pi[inversion] / pi[genome])) +
    labs(y = ylab) + theme_nature() + theme(legend.position='none', axis.text.x = element_text(size=7), axis.text.y = element_text(size=7))
  if (!is.null(hline1)) p <- p + geom_hline(yintercept = hline1, linetype = 'dashed', linewidth = 0.35)
  p
}

summarise_by_tebin <- function(df, value_col, te_breaks = seq(0,1,by=0.2)) {
  v <- enquo(value_col); te_mids <- (te_breaks[-length(te_breaks)] + te_breaks[-1]) / 2
  dat <- df %>% filter(!is.na(te_prop_union), te_prop_union >= 0, te_prop_union <= 1, !is.na(!!v), is.finite(!!v)) %>% mutate(te_mid = te_mids[cut(te_prop_union, breaks = te_breaks, labels = FALSE, include.lowest = TRUE, right = FALSE)]) %>% filter(!is.na(te_mid))
  summary_tab <- summarise_boot_ci(dat, !!v, group_vars = c('kingdom','te_mid'))
  stat_tab <- dat %>% group_by(te_mid) %>% filter(n_distinct(kingdom)==2) %>% wilcox_test(as.formula(paste0(rlang::as_name(v),' ~ kingdom'))) %>% adjust_pvalue(method='BH') %>% add_significance('p.adj') %>% ungroup()
  bin_y <- summary_tab %>% group_by(te_mid) %>% summarise(y = max(ci_high, na.rm=TRUE), .groups='drop') %>% mutate(y = y * 1.06)
  list(summary = summary_tab, stats = left_join(stat_tab, bin_y, by='te_mid'))
}

plot_tebin_trait <- function(sum_stat_list, ylab, hline1 = 1) {
  summary_tab <- sum_stat_list$summary; stat_tab <- sum_stat_list$stats
  p <- ggplot(summary_tab, aes(x = te_mid, y = mean, colour = kingdom, group = kingdom)) +
    geom_line(linewidth = 0.7) + geom_point(size = 1.4) + geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.03, linewidth = 0.4) +
    geom_text(data = stat_tab, aes(x = te_mid, y = y, label = p.adj.signif), inherit.aes = FALSE, size = 2.3) +
    scale_color_manual(values = kingdom_cols) + scale_x_continuous(breaks = seq(0.1,0.9,0.2), labels = scales::percent(seq(0.1,0.9,0.2), accuracy = 1), name = 'TE proportion in inversion') +
    labs(y = ylab) + theme_nature() + theme(legend.position='none', axis.text.x = element_text(size=7), axis.text.y = element_text(size=7))
  if (!is.null(hline1)) p <- p + geom_hline(yintercept = hline1, linetype = 'dashed', linewidth = 0.35)
  p
}

plot_lenclass_box <- function(df, value_col, ylab, hline1 = NULL, ylim_q = c(0.01,0.99), dodge_w = 0.65, jitter_w = 0.15, alpha_j = 0.12) {
  v <- enquo(value_col); dat <- df %>% filter(!is.na(size_class), !is.na(!!v), is.finite(!!v))
  yq <- quantile(pull(dat, !!v), ylim_q, na.rm = TRUE)
  stat_tab <- dat %>% group_by(size_class) %>% filter(n_distinct(kingdom)==2) %>% wilcox_test(as.formula(paste0(rlang::as_name(v),' ~ kingdom'))) %>% adjust_pvalue(method='BH') %>% add_significance('p.adj') %>% ungroup() %>% mutate(y = as.numeric(yq[2]) * 1.06)
  p <- ggplot(dat, aes(x = size_class, y = !!v, fill = kingdom)) +
    geom_jitter(aes(color = kingdom), position = position_jitterdodge(jitter.width = jitter_w, dodge.width = dodge_w), size = 0.4, alpha = alpha_j, show.legend = FALSE) +
    geom_boxplot(position = position_dodge(width = dodge_w), width = 0.45, linewidth = 0.4, outlier.shape = NA, alpha = 0.85) +
    scale_fill_manual(values = kingdom_cols) + scale_color_manual(values = kingdom_cols) + coord_cartesian(ylim = c(as.numeric(yq[1]), as.numeric(yq[2]) * 1.12)) +
    labs(x = 'Inversion length class', y = ylab) + theme_nature() + theme(legend.position='none', axis.text.x = element_text(size=7)) +
    geom_text(data = stat_tab, aes(x = size_class, y = y, label = p.adj.signif), inherit.aes = FALSE, size = 2.3)
  if (!is.null(hline1)) p <- p + geom_hline(yintercept = hline1, linetype = 'dashed', linewidth = 0.35)
  p
}

plot_density_kingdom_ratio <- function(df, xcol, xlab, xlim = c(0,2)) {
  dat <- df %>% filter(is.finite(.data[[xcol]]), !is.na(.data[[xcol]]))
  wt <- tryCatch(wilcox.test(dat[[xcol]] ~ dat$kingdom, exact = FALSE)$p.value, error = function(e) NA_real_)
  lab <- ifelse(is.na(wt), 'p = NA', ifelse(wt < 0.001, 'p < 0.001***', paste0('p = ', signif(wt,2), ' ', p_to_star(wt))))
  dens_max <- safe_density_max(dat[[xcol]])
  med_tbl <- dat %>% group_by(kingdom) %>% summarise(med = median(.data[[xcol]], na.rm = TRUE), .groups='drop')
  ggplot(dat, aes(x = .data[[xcol]], color = kingdom, fill = kingdom)) +
    geom_density(alpha = 0.25, adjust = 1.0, linewidth = 0.7) +
    geom_vline(data = med_tbl, aes(xintercept = med, color = kingdom), linewidth = 0.6, linetype = 'dashed', show.legend = FALSE) +
    scale_color_manual(values = kingdom_cols) + scale_fill_manual(values = kingdom_cols) + coord_cartesian(xlim = xlim) +
    annotate('text', x = xlim[2] - 0.02 * diff(xlim), y = dens_max * 0.75, label = lab, hjust = 1, vjust = 1, size = 2.6, color = 'black') +
    labs(x = xlab, y = 'Density') + theme_nature() + theme(legend.position = c(0.98,0.98), legend.justification = c(1,1), legend.direction='vertical', legend.text = element_text(size=6), legend.key.size = grid::unit(3,'mm'))
}

plot_density_size <- function(df, xcol, xlab, xlim = c(0,2), facet_by_kingdom = TRUE) {
  dat <- df %>% filter(is.finite(.data[[xcol]]), !is.na(size_bin))
  kw <- tryCatch(kruskal.test(dat[[xcol]] ~ dat$size_bin)$p.value, error = function(e) NA_real_)
  lab <- ifelse(is.na(kw), 'p = NA', ifelse(kw < 0.001, 'p < 0.001***', paste0('p = ', signif(kw,2), ' ', p_to_star(kw))))
  dens_max <- safe_density_max(dat[[xcol]])
  med_size <- dat %>% group_by(size_bin) %>% summarise(med = median(.data[[xcol]], na.rm = TRUE), .groups='drop')
  size_cols <- c('small'='#9ecae1','medium'='#fdae6b','large'='#a1d99b')
  p <- ggplot(dat, aes(x = .data[[xcol]], color = size_bin, fill = size_bin)) +
    geom_density(alpha = 0.20, adjust = 1.0, linewidth = 0.7) +
    geom_vline(data = med_size, aes(xintercept = med, color = size_bin), linewidth = 0.6, linetype = 'dashed', show.legend = FALSE) +
    scale_color_manual(values = size_cols) + scale_fill_manual(values = size_cols) + coord_cartesian(xlim = xlim) +
    annotate('text', x = xlim[2] - 0.02 * diff(xlim), y = dens_max * 0.95, label = lab, hjust = 1, vjust = 1, size = 2.6, color = 'black') +
    labs(x = xlab, y = 'Density') + theme_nature() + theme(legend.position = c(0.98,0.98), legend.justification = c(1,1), legend.direction='vertical', legend.text = element_text(size=6), legend.key.size = grid::unit(3,'mm'))
  if (facet_by_kingdom) p <- p + facet_wrap(~ kingdom)
  p
}
