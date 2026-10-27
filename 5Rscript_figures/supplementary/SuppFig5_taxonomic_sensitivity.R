#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(scales)
})

set.seed(20260925)
get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_flag <- '--file='
  path <- sub(file_flag, '', args[grepl(file_flag, args)])
  if (length(path) > 0) dirname(normalizePath(path[1])) else getwd()
}

script_dir <- get_script_dir()
PROJ_ROOT <- if (basename(script_dir) %in% c('extended_data','supplementary','main_figures')) dirname(script_dir) else script_dir

OUT_DIR <- file.path(PROJ_ROOT, 'outputs', 'supplementary')
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

PER_GENOME <- "per_genome.csv"
PER_INV <- "per_inversion.csv"

pg <- read_csv(PER_GENOME, show_col_types = FALSE, progress = FALSE)
pi <- read_csv(PER_INV, show_col_types = FALSE, progress = FALSE)

required_pg <- c(
  "sample_id", "Species", "kingdom", "order", "major_clade", "gene_counts",
  "inv_count_total", "inv_count_small_lt10kb", "inv_count_medium_10kb_1Mb",
  "inv_count_large_gt1Mb", "inv_count_per_Gb", "inv_len_mean_bp",
  "inv_len_median_bp", "TE_total_pct", "gene_density", "asm_total_len_bp",
  "pi_genome_wide"
)
required_pi <- c("inv_id", "sample_id", "kingdom", "inv_len_bp")
stopifnot(all(required_pg %in% names(pg)), all(required_pi %in% names(pi)))

vertebrate_clades <- c("Birds", "Fishes", "Mammals", "Reptiles")

group_cols <- c(
  "Plants" = "#009E73",
  "Insects" = "#E69F00",
  "Vertebrates" = "#0072B2",
  "Other non-insect Metazoa" = "#999999"
)
comparison_cols <- c(
  "All Metazoa vs Plants" = "#222222",
  "Non-insect Metazoa vs Plants" = "#777777",
  "Vertebrates vs Plants" = "#0072B2"
)
x_labels <- c(
  "Plants" = "Plants",
  "Insects" = "Insects",
  "Vertebrates" = "Vertebrates",
  "Other non-insect Metazoa" = "Other non-insect\nMetazoa"
)

pg2 <- pg %>%
  mutate(
    inv_per_10k_genes = if_else(gene_counts > 0 & is.finite(gene_counts) & is.finite(inv_count_total), inv_count_total / gene_counts * 10000, NA_real_),
    log10_total_count = if_else(inv_count_total > 0 & is.finite(inv_count_total), log10(inv_count_total), NA_real_),
    log10_inv_per_10k_genes = if_else(inv_per_10k_genes > 0 & is.finite(inv_per_10k_genes), log10(inv_per_10k_genes), NA_real_),
    log10_inv_count_per_Gb = if_else(inv_count_per_Gb > 0 & is.finite(inv_count_per_Gb), log10(inv_count_per_Gb), NA_real_),
    log10_median_inv_length = if_else(inv_len_median_bp > 0 & is.finite(inv_len_median_bp), log10(inv_len_median_bp), NA_real_),
    log10_mean_inv_length = if_else(inv_len_mean_bp > 0 & is.finite(inv_len_mean_bp), log10(inv_len_mean_bp), NA_real_),
    prop_large = if_else(inv_count_total > 0 & is.finite(inv_count_total) & is.finite(inv_count_large_gt1Mb), inv_count_large_gt1Mb / inv_count_total, NA_real_),
    sensitivity_group = case_when(
      kingdom == "Viridiplantae" ~ "Plants",
      kingdom == "Metazoa" & major_clade == "Insects" ~ "Insects",
      kingdom == "Metazoa" & major_clade %in% vertebrate_clades ~ "Vertebrates",
      kingdom == "Metazoa" & !major_clade %in% c("Insects", vertebrate_clades) ~ "Other non-insect Metazoa",
      TRUE ~ NA_character_
    ),
    sensitivity_group = factor(sensitivity_group, levels = c("Plants", "Insects", "Vertebrates", "Other non-insect Metazoa"))
  )

pi2 <- pi %>%
  select(inv_id, sample_id, kingdom, inv_len_bp) %>%
  inner_join(pg2 %>% select(sample_id, Species, major_clade, sensitivity_group), by = "sample_id") %>%
  filter(!is.na(sensitivity_group), is.finite(inv_len_bp), inv_len_bp >= 200) %>%
  mutate(log10_inv_len = log10(inv_len_bp))

theme_nature <-theme_bw(base_size = 8.0) +
    theme(
      legend.title = element_blank(),
      panel.grid      = element_blank(),
      axis.text       = element_text(color = "black"),
      axis.title      = element_text(color = "black"),
      axis.ticks      = element_line(color = "black"),
      panel.border    = element_rect(color = "black", linewidth = 0.5),
      legend.position = "none",
      axis.text.x = element_text(angle = 24, hjust = 1, vjust = 1),
      plot.title      = element_blank(),
      plot.tag        = element_text(face = "bold", size = 8)
    )

compact_count_labels <- function(x) {
  v <- 10^x
  case_when(
    v >= 1e6 ~ paste0(signif(v / 1e6, 3), "M"),
    v >= 1e3 ~ paste0(signif(v / 1e3, 3), "k"),
    v >= 1 ~ format(signif(v, 3), trim = TRUE, scientific = FALSE),
    TRUE ~ format(signif(v, 2), trim = TRUE, scientific = TRUE)
  )
}

compact_bp_labels <- function(x) {
  v <- 10^x
  case_when(
    v >= 1e6 ~ paste0(signif(v / 1e6, 3), " Mb"),
    v >= 1e3 ~ paste0(signif(v / 1e3, 3), " kb"),
    TRUE ~ paste0(signif(v, 3), " bp")
  )
}

compact_log_breaks <- function(x, max_labels = 3L) {
  x <- x[is.finite(x)]
  if (!length(x)) return(numeric())
  lo <- floor(min(x))
  hi <- ceiling(max(x))
  ints <- seq(lo, hi)
  if (length(ints) <= max_labels) return(ints)
  step <- ceiling((length(ints) - 1) / (max_labels - 1))
  br <- seq(lo, hi, by = step)
  if (tail(br, 1) != hi && length(br) < max_labels) br <- c(br, hi)
  br
}

backtransform_log_axis <- function(p, values, unit = c("count", "bp"), max_labels = 3L) {
  unit <- match.arg(unit)
  br <- compact_log_breaks(values, max_labels = max_labels)
  lab_fun <- if (unit == "bp") compact_bp_labels else compact_count_labels
  p + scale_y_continuous(breaks = br, labels = lab_fun(br))
}
make_boxplot <- function(data, yvar, ylab, percent_axis = FALSE, natural_unit = NULL) {
  p <- data %>%
    filter(!is.na(sensitivity_group), is.finite(.data[[yvar]])) %>%
    ggplot(aes(sensitivity_group, .data[[yvar]], fill = sensitivity_group)) +
    geom_boxplot(width = 0.58, outlier.shape = NA, linewidth = 0.38, alpha = 0.72) +
    geom_jitter(aes(colour = sensitivity_group), width = 0.15, size = 0.68, alpha = 0.28) +
    scale_fill_manual(values = group_cols) +
    scale_colour_manual(values = group_cols) +
    scale_x_discrete(labels = x_labels) +
    labs(x = NULL, y = ylab) +
    guides(fill = "none", colour = "none") +
    theme_nature
  if (percent_axis) p <- p + scale_y_continuous(labels = percent_format(accuracy = 1))
  if (!is.null(natural_unit)) p <- backtransform_log_axis(p, data[[yvar]], unit = natural_unit, max_labels = 4L)
  p
}

pA <- make_boxplot(pg2, "log10_total_count", "Total inversion count", natural_unit="count")
pB <- make_boxplot(pg2, "log10_inv_per_10k_genes", "Inversions per 10,000 genes", natural_unit="count")
pC <- make_boxplot(pg2, "log10_inv_count_per_Gb", "Inversion count per Gb", natural_unit="count")
pD <- make_boxplot(pg2, "log10_median_inv_length", "Median inversion length", natural_unit="bp")
pE <- make_boxplot(pg2, "log10_mean_inv_length", "Mean inversion length", natural_unit="bp")
pF <- make_boxplot(pg2, "prop_large", "Proportion of inversions >1 Mb", TRUE)

pG <- pi2 %>%
  ggplot(aes(log10_inv_len, colour = sensitivity_group, fill = sensitivity_group)) +
  geom_density(linewidth = 0.72, alpha = 0.09) +
  scale_colour_manual(values = group_cols) +
  scale_fill_manual(values = group_cols) +
  scale_x_continuous(breaks=c(2,4,6,8), labels=c("100 bp","10 kb","1 Mb","100 Mb")) +
  labs(x = "Inversion length", y = "Density") +
  theme_nature +
  theme(
    legend.background = element_rect(fill = "transparent", colour = NA),
    axis.text.y = element_blank(), axis.ticks.y = element_blank(),
    legend.position = c(0.70, 0.85), legend.text = element_text(size = 5.5),
    legend.key.size = unit(3, "mm")
  )

effect_metrics <- tribble(
  ~metric, ~metric_label,
  "inv_count_total", "Total inversion count",
  "inv_per_10k_genes", "Inversions / 10k genes",
  "inv_count_per_Gb", "Inversions / Gb",
  "inv_len_median_bp", "Median inversion length",
  "inv_len_mean_bp", "Mean inversion length",
  "prop_large", "Proportion >1 Mb"
)

bootstrap_log2_median_ratio <- function(animal_values, plant_values, B = 10000, seed = 1) {
  a <- animal_values[is.finite(animal_values)]
  p <- plant_values[is.finite(plant_values)]
  ma <- median(a)
  mp <- median(p)
  if (!is.finite(ma) || !is.finite(mp) || ma <= 0 || mp <= 0) stop("Group median <= 0; cannot calculate log2 median ratio.")
  set.seed(seed)
  boot <- replicate(B, {
    ba <- median(sample(a, length(a), replace = TRUE))
    bp <- median(sample(p, length(p), replace = TRUE))
    if (ba > 0 && bp > 0) log2(ba / bp) else NA_real_
  })
  q <- quantile(boot, c(0.025, 0.975), na.rm = TRUE, names = FALSE)
  tibble(
    animal_median = ma,
    plant_median = mp,
    median_ratio = ma / mp,
    ratio_ci_low = 2^q[1],
    ratio_ci_high = 2^q[2],
    log2_ratio = log2(ma / mp),
    ci_low = q[1],
    ci_high = q[2]
  )
}

plants <- pg2 %>% filter(kingdom == "Viridiplantae")
comparison_sets <- list(
  "All Metazoa vs Plants" = pg2 %>% filter(kingdom == "Metazoa"),
  "Non-insect Metazoa vs Plants" = pg2 %>% filter(kingdom == "Metazoa", !(major_clade %in% "Insects")),
  "Vertebrates vs Plants" = pg2 %>% filter(kingdom == "Metazoa", major_clade %in% vertebrate_clades)
)

effect_results <- imap_dfr(comparison_sets, function(animal_df, comp_name) {
  map_dfr(seq_len(nrow(effect_metrics)), function(i) {
    metric <- effect_metrics$metric[i]
    bootstrap_log2_median_ratio(animal_df[[metric]], plants[[metric]], seed = 20260925 + i + match(comp_name, names(comparison_sets)) * 100) %>%
      mutate(
        comparison = comp_name,
        metric = metric,
        metric_label = effect_metrics$metric_label[i],
        n_animal = sum(is.finite(animal_df[[metric]])),
        n_plant = sum(is.finite(plants[[metric]]))
      )
  })
}) %>%
  mutate(
    comparison = factor(comparison, levels = names(comparison_sets)),
    metric_label = factor(metric_label, levels = rev(effect_metrics$metric_label))
  )

pH <- effect_results %>%
  ggplot(aes(metric_label, log2_ratio, colour = comparison, group = comparison)) +
  geom_hline(yintercept = 0, linetype = 2, linewidth = 0.35, colour = "grey45") +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.15, linewidth = 0.48, position = position_dodge(0.58)) +
  geom_point(size = 1.65, position = position_dodge(0.58)) +
  scale_colour_manual(
    values = comparison_cols,
    labels = c("All Metazoa", "Non-insect Metazoa", "Vertebrates")
  ) +
  coord_flip() +
  labs(x = NULL, y = expression(log[2] * "(animal / plant median)")) +
  theme_nature +
  theme(legend.position = c(0.75, 0.15), 
        legend.background = element_rect(fill = "transparent", colour = NA),
        legend.text = element_text(size = 5.5),
        legend.key.size = unit(3, "mm"),
        axis.text.y = element_text(size = 6.8))

insects <- pg2 %>% filter(kingdom == "Metazoa", major_clade == "Insects")
other_metazoa <- pg2 %>% filter(kingdom == "Metazoa", !(major_clade %in% "Insects"))
stopifnot(nrow(insects) >= 50)

set.seed(20260925)
downsample_results <- map_dfr(seq_len(1000), function(rep_id) {
  animal_rep <- bind_rows(slice_sample(insects, n = 50), other_metazoa)
  map_dfr(seq_len(nrow(effect_metrics)), function(i) {
    metric <- effect_metrics$metric[i]
    a <- animal_rep[[metric]][is.finite(animal_rep[[metric]])]
    p <- plants[[metric]][is.finite(plants[[metric]])]
    ma <- median(a)
    mp <- median(p)
    tibble(
      replicate = rep_id,
      metric = metric,
      metric_label = effect_metrics$metric_label[i],
      log2_ratio = if (ma > 0 && mp > 0) log2(ma / mp) else NA_real_
    )
  })
})

downsample_summary <- downsample_results %>%
  group_by(metric, metric_label) %>%
  summarise(
    median_log2_ratio = median(log2_ratio, na.rm = TRUE),
    interval_low = quantile(log2_ratio, 0.025, na.rm = TRUE),
    interval_high = quantile(log2_ratio, 0.975, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(metric_label = factor(metric_label, levels = rev(effect_metrics$metric_label)))

pI <- downsample_summary %>%
  ggplot(aes(metric_label, median_log2_ratio)) +
  geom_hline(yintercept = 0, linetype = 2, linewidth = 0.35, colour = "grey45") +
  geom_errorbar(aes(ymin = interval_low, ymax = interval_high), width = 0.16, linewidth = 0.52) +
  geom_point(size = 1.8) +
  coord_flip() +
  labs(x = NULL, y = expression(log[2] * "(animal / plant median)")) +
  theme_nature +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
        axis.text.y = element_blank() )


final_fig <- ((pA | pB | pC) / (pD | pE | pF) / (pG | pH | pI)) +
  plot_layout(heights = c(1, 1, 1)) +
  plot_annotation(tag_levels = "a")

ggsave(file.path(OUT_DIR, "SuppFig5_taxonomic_sensitivity.pdf"), final_fig, width = 210, height = 220, units = "mm", device = cairo_pdf)
ggsave(file.path(OUT_DIR, "SuppFig5_taxonomic_sensitivity.png"), final_fig, width = 210, height = 220, units = "mm", dpi = 300, bg = "white")

response_table <- effect_results %>%
  select(comparison, metric_label, n_animal, n_plant, animal_median, plant_median, median_ratio, ratio_ci_low, ratio_ci_high, log2_ratio, ci_low, ci_high) %>%
  arrange(metric_label, comparison)

write_csv(response_table, file.path(OUT_DIR, "SuppFig5_response_letter_effect_sizes_ALL_METRICS.csv"))
write_csv(downsample_summary, file.path(OUT_DIR, "SuppFig5_insect_downsampling_summary.csv"))
