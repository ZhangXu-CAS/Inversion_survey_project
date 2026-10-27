############################################################
## Extended Data Figure 3
## Per-genome distributions across major clades
##
## Panel order:
##   a  Inversion count per Gb
##   b  Total inversion count per genome
##   c  Median inversion length
##   d  Proportion of large inversions (>1 Mb)
############################################################


get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_flag <- '--file='
  path <- sub(file_flag, '', args[grepl(file_flag, args)])
  if (length(path) > 0) dirname(normalizePath(path[1])) else getwd()
}
script_dir <- get_script_dir()
PROJ_ROOT <- if (basename(script_dir) %in% c('extended_data','supplementary','main_figures')) dirname(script_dir) else script_dir
source(file.path(PROJ_ROOT, 'helpers', 'figure_plot_helpers.R'))
OUT_DIR <- file.path(PROJ_ROOT, 'outputs', 'extended_data')
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
IN_PER_GENOME <- file.path(PROJ_ROOT, 'per_genome.csv')
IN_PER_INVERSION <- file.path(PROJ_ROOT, 'per_inversion.csv')

LARGE_INV_BP <- 1e6

per_genome <- load_per_genome() %>%
  filter(kingdom %in% KINGDOM_LEVELS_MV) %>%
  droplevels()

require_columns(
  per_genome,
  c("major_clade", "sample_id", "inv_count_total", "inv_count_per_Gb", "median_log10_inv_length"),
  "per_genome"
)

plot_panel_clade <- function(df, y, ylab, y_is_prop = FALSE, natural_unit = NULL) {
  y <- enquo(y)
  dat <- df %>% filter(!is.na(!!y), is.finite(!!y))

  p <- ggplot(dat, aes(x = clade_key, y = !!y)) +
    geom_boxplot(
      aes(fill = kingdom),
      width = 0.55,
      linewidth = 0.35,
      outlier.shape = NA,
      alpha = 0.85
    ) +
    geom_jitter(
      aes(color = kingdom),
      width = 0.15,
      size = 0.35,
      alpha = 0.15,
      show.legend = FALSE
    ) +
    facet_wrap(~ kingdom, scales = "free_x", nrow = 1) +
    scale_fill_manual(values = kingdom_cols) +
    scale_color_manual(values = kingdom_cols) +
    scale_x_discrete(labels = function(x) sub("^.*?:", "", x)) +
    labs(x = NULL, y = ylab) +
    theme_nature() +
    theme(
      axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1, size = 5),
      axis.title.y = element_text(size = 6)
    )

  if (y_is_prop) p <- p + scale_y_continuous(labels = percent_format(accuracy = 1))
  if (!is.null(natural_unit)) p <- backtransform_log_axis(p, pull(dat, !!y), axis="y", unit=natural_unit, max_labels=4)
  p
}

## derive proportion of large inversions (>1 Mb)
ge1mb_col <- find_ge1mb_col(names(per_genome))

if (!is.na(ge1mb_col)) {
  per_genome <- per_genome %>%
    mutate(
      inv_count_ge_1Mb = as.numeric(.data[[ge1mb_col]]),
      prop_large_inv_1Mb = if_else(inv_count_total > 0, inv_count_ge_1Mb / inv_count_total, NA_real_)
    )
} else {
  per_inv <- load_per_inversion()
  require_columns(per_inv, c("sample_id", "inv_len_bp"), "per_inversion")

  pinv_sum <- per_inv %>%
    mutate(inv_len_bp = as.numeric(inv_len_bp)) %>%
    filter(!is.na(inv_len_bp), inv_len_bp > 0) %>%
    group_by(sample_id) %>%
    summarise(inv_count_ge_1Mb = sum(inv_len_bp > LARGE_INV_BP), .groups = "drop")

  per_genome <- per_genome %>%
    left_join(pinv_sum, by = "sample_id") %>%
    mutate(
      inv_count_ge_1Mb = replace_na(inv_count_ge_1Mb, 0),
      prop_large_inv_1Mb = if_else(inv_count_total > 0, inv_count_ge_1Mb / inv_count_total, NA_real_)
    )
}

per_genome <- per_genome %>%
  mutate(
    kingdom = factor(kingdom, levels = KINGDOM_LEVELS_MV),
    clade_key = paste0(as.character(kingdom), ":", major_clade),
    log10_inv_count_per_Gb = if_else(is.finite(inv_count_per_Gb) & inv_count_per_Gb > 0,
                                     safe_log10(inv_count_per_Gb), NA_real_),
    log10_inv_count_total = if_else(is.finite(inv_count_total) & inv_count_total > 0,
                                    safe_log10(inv_count_total), NA_real_)
  )

ed3a <- plot_panel_clade(
  per_genome,
  log10_inv_count_per_Gb,
  ylab = "Inversion count per Gb",
  natural_unit = "count"
)

ed3b <- plot_panel_clade(
  per_genome,
  log10_inv_count_total,
  ylab = "Total inversion count per genome",
  natural_unit = "count"
)

ed3c <- plot_panel_clade(
  per_genome,
  median_log10_inv_length,
  ylab = "Median inversion length",
  natural_unit = "bp"
)

ed3d <- plot_panel_clade(
  per_genome,
  prop_large_inv_1Mb,
  ylab = "Proportion of large inversions",
  y_is_prop = TRUE
)

Extended_Data_Fig3 <- (ed3a | ed3b) / (ed3c | ed3d) +
  plot_annotation(tag_levels = "a")

save_pdf(Extended_Data_Fig3, "Extended_Data_Fig3.pdf", 180, 120)
save_png(Extended_Data_Fig3, "Extended_Data_Fig3.png", 180, 120, dpi = 300)

message("Done: Extended_Data_Fig3.pdf / Extended_Data_Fig3.png")
