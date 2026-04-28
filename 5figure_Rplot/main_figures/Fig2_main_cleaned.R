############################################################
## Figure 2 (main)
## Inversion number and length across kingdoms
##
## Panel order (matches manuscript):
##   a  Inversion count per Gb
##   b  Total inversion count per genome
##   c  Genome-wide π vs inversion count per Gb
##   d  Median inversion length per genome
##   e  Mean inversion length per genome
##   f  Genome-wide π vs median inversion length
##   g  Per-inversion log10 length density
##   h  Mean size-class composition (small / medium / large)
############################################################


get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_flag <- '--file='
  path <- sub(file_flag, '', args[grepl(file_flag, args)])
  if (length(path) > 0) dirname(normalizePath(path[1])) else getwd()
}
SCRIPT_DIR <- get_script_dir()
PROJ_ROOT <- normalizePath(file.path(SCRIPT_DIR, '..'), winslash = '/', mustWork = FALSE)
source(file.path(PROJ_ROOT, 'helpers', 'figure_plot_helpers.R'))
OUT_DIR <- file.path(PROJ_ROOT, 'outputs', 'main_figures')
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

per_genome <- load_per_genome()
per_inv <- load_per_inversion()

require_columns(
  per_genome,
  c("log10_inv_count_per_Gb", "log10_inv_count_total", "median_log10_inv_length", "mean_log10_inv_length", "pi_genome_wide"),
  "per_genome"
)

if (!all(c("prop_small", "prop_medium", "prop_large") %in% names(per_genome))) {
  stop("Figure 2 panel h requires size-class counts in per_genome.csv to derive prop_small / prop_medium / prop_large.", call. = FALSE)
}

## a
fig2a_dat <- per_genome %>% filter(!is.na(log10_inv_count_per_Gb))
fig2a <- boxpanel_kingdom(
  fig2a_dat,
  log10_inv_count_per_Gb,
  ylab = expression(paste("Inversion count per Gb (log"[10], ")")),
  x_labels = make_x_labels_with_n(fig2a_dat)
)

## b
fig2b_dat <- per_genome %>% filter(!is.na(log10_inv_count_total))
fig2b <- boxpanel_kingdom(
  fig2b_dat,
  log10_inv_count_total,
  ylab = expression(paste("Total inversion count (log"[10], ")")),
  x_labels = make_x_labels_with_n(fig2b_dat)
)

## c
fig2c <- corr_panel_kingdom(
  per_genome,
  xvar = "pi_genome_wide",
  yvar = "inv_count_per_Gb",
  xlab = "Genome-wide nucleotide diversity (π)",
  ylab = expression(paste("Inversion count per Gb (log"[10], ")")),
  logy = TRUE,
  method_stats = "pearson",
  legend_pos = c(0.70, 0.20),
  legend_text_size = 5
)

## d
fig2d_dat <- per_genome %>% filter(!is.na(median_log10_inv_length))
fig2d <- boxpanel_kingdom(
  fig2d_dat,
  median_log10_inv_length,
  ylab = expression(paste("Median inversion length (log"[10], " bp)")),
  x_labels = make_x_labels_with_n(fig2d_dat)
)

## e
fig2e_dat <- per_genome %>% filter(!is.na(mean_log10_inv_length))
fig2e <- boxpanel_kingdom(
  fig2e_dat,
  mean_log10_inv_length,
  ylab = expression(paste("Mean inversion length (log"[10], " bp)")),
  x_labels = make_x_labels_with_n(fig2e_dat)
)

## f
fig2f <- corr_panel_kingdom(
  per_genome,
  xvar = "pi_genome_wide",
  yvar = "median_log10_inv_length",
  xlab = "Genome-wide nucleotide diversity (π)",
  ylab = expression(paste("Median inversion length (log"[10], " bp)")),
  logy = FALSE,
  method_stats = "pearson",
  legend_pos = c(0.70, 0.80),
  legend_text_size = 5
)

## g
fig2g <- ggplot(per_inv, aes(x = log10_len, color = kingdom, fill = kingdom)) +
  geom_density(alpha = 0.25, adjust = 1.0) +
  scale_color_manual(values = kingdom_cols) +
  scale_fill_manual(values = kingdom_cols) +
  labs(x = expression(paste("Inversion length (log"[10], " bp)")), y = "Density") +
  theme_nature() +
  theme(
    legend.position = c(0.98, 0.98),
    legend.justification = c(1, 1),
    legend.direction = "vertical",
    legend.title = element_blank(),
    legend.text = element_text(size = 6),
    legend.key.size = grid::unit(3, "mm"),
    plot.margin = margin(1.5, 3, 1.5, 2, unit = "mm")
  )

## h
size_summary <- per_genome %>%
  filter(!is.na(prop_small), !is.na(prop_medium), !is.na(prop_large)) %>%
  pivot_longer(
    cols = c(prop_small, prop_medium, prop_large),
    names_to = "size_class",
    values_to = "prop"
  ) %>%
  mutate(
    size_class = recode(
      size_class,
      prop_small = "small",
      prop_medium = "medium",
      prop_large = "large"
    ),
    size_class = factor(size_class, levels = c("small", "medium", "large"))
  ) %>%
  group_by(kingdom, size_class) %>%
  summarise(mean_prop = mean(prop, na.rm = TRUE), .groups = "drop")

fig2h <- ggplot(size_summary, aes(x = kingdom, y = mean_prop, fill = size_class)) +
  geom_col(width = 0.7, color = "black", linewidth = 0.2) +
  scale_fill_manual(values = size_cols) +
  scale_x_discrete(labels = make_x_labels_plain(levels(per_genome$kingdom))) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(x = NULL, y = "Mean proportion of inversions") +
  theme_nature() +
  theme(
    legend.position = "right",
    legend.direction = "vertical",
    legend.title = element_blank(),
    legend.text = element_text(size = 6),
    legend.key.size = grid::unit(3, "mm"),
    plot.margin = margin(1.5, 2, 1.5, 2, unit = "mm")
  )

Figure2_main <- wrap_plots(
  fig2a, fig2b, fig2c,
  fig2d, fig2e, fig2f,
  fig2g, fig2h,
  ncol = 3
) +
  plot_annotation(tag_levels = "a")

save_pdf(Figure2_main, "Fig2_main.pdf", 180, 170)
save_png(Figure2_main, "Fig2_main.png", 180, 170, dpi = 300)

message("Done: Fig2_main.pdf / Fig2_main.png")
