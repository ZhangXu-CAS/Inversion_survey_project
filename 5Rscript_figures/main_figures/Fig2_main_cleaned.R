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
script_dir <- get_script_dir()
PROJ_ROOT <- if (basename(script_dir) == 'main_figures') dirname(script_dir) else script_dir

source(file.path(PROJ_ROOT, 'helpers', 'figure_plot_helpers.R'))
OUT_DIR <- file.path(PROJ_ROOT, 'outputs', 'main_figures')
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
IN_PER_GENOME <- file.path(PROJ_ROOT, 'per_genome.csv')
IN_PER_INVERSION <- file.path(PROJ_ROOT, 'per_inversion.csv')

per_genome <- load_per_genome()
per_inv <- load_per_inversion()
PGLS_FILE <- file.path(PROJ_ROOT, 'Fig2_PGLS_beta_95CI.csv')
pgls_c_labels <- read_pgls_labels(PGLS_FILE, 'log10_inv_count_per_Gb')
pgls_f_labels <- read_pgls_labels(PGLS_FILE, 'median_log10_inv_length')

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
  ylab = "Inversion count per Gb",
  x_labels = make_x_labels_with_n(fig2a_dat),
  natural_unit = 'count'
)

## b
fig2b_dat <- per_genome %>% filter(!is.na(log10_inv_count_total))
fig2b <- boxpanel_kingdom(
  fig2b_dat,
  log10_inv_count_total,
  ylab = "Total inversion count",
  x_labels = make_x_labels_with_n(fig2b_dat),
  natural_unit = 'count'
)

## c
fig2c <- corr_panel_kingdom(
  per_genome,
  xvar = "pi_genome_wide",
  yvar = "inv_count_per_Gb",
  xlab = "Genome-wide haplotypic divergence (π)",
  ylab = "Inversion count per Gb",
  logy = TRUE,
  method_stats = "pearson",
  legend_pos = c(0.65, 0.20),
  legend_text_size = 5,
  legend_labels_override = pgls_c_labels,
  natural_y_unit = 'count'
)

## d
fig2d_dat <- per_genome %>% filter(!is.na(median_log10_inv_length))
fig2d <- boxpanel_kingdom(
  fig2d_dat,
  median_log10_inv_length,
  ylab = "Median inversion length",
  x_labels = make_x_labels_with_n(fig2d_dat),
  natural_unit = 'bp'
)

## e
fig2e_dat <- per_genome %>% filter(!is.na(mean_log10_inv_length))
fig2e <- boxpanel_kingdom(
  fig2e_dat,
  mean_log10_inv_length,
  ylab = "Mean inversion length",
  x_labels = make_x_labels_with_n(fig2e_dat),
  natural_unit = 'bp'
)

## f
fig2f <- corr_panel_kingdom(
  per_genome,
  xvar = "pi_genome_wide",
  yvar = "median_log10_inv_length",
  xlab = "Genome-wide haplotypic divergence (π)",
  ylab = "Median inversion length",
  logy = FALSE,
  method_stats = "pearson",
  legend_pos = c(0.65, 0.80),
  legend_text_size = 5,
  legend_labels_override = pgls_f_labels,
  natural_y_unit = 'bp'
)

## g
fig2g <- ggplot(per_inv, aes(x = log10_len, color = kingdom, fill = kingdom)) +
  geom_density(alpha = 0.25, adjust = 1.0) +
  scale_color_manual(values = kingdom_cols) +
  scale_fill_manual(values = kingdom_cols) +
  scale_x_continuous(breaks = c(2,4,6,8), labels = c('100 bp','10 kb','1 Mb','100 Mb')) +
  labs(x = "Inversion length", y = "Density") +
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
      prop_small = "Small",
      prop_medium = "Medium",
      prop_large = "Large"
    ),
    size_class = factor(size_class, levels = c("Small", "Medium", "Large"))
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


row3 <- (fig2g | fig2h | plot_spacer()) +
plot_layout(widths = c(1.5, 1, 1e-6), guides = "keep")

Figure2_main <- (( fig2a |fig2b | fig2c) /
           ( fig2d |fig2e | fig2f)/
           row3 ) +
  plot_annotation(tag_levels = "a")

fig2_effects <- bind_rows(
  bootstrap_diff_ci(fig2a_dat, 'log10_inv_count_per_Gb') %>% mutate(panel='a', metric='log10 inversion count per Gb'),
  bootstrap_diff_ci(fig2b_dat, 'log10_inv_count_total') %>% mutate(panel='b', metric='log10 total inversion count'),
  bootstrap_diff_ci(fig2d_dat, 'median_log10_inv_length') %>% mutate(panel='d', metric='genome-median log10 inversion length'),
  bootstrap_diff_ci(fig2e_dat, 'mean_log10_inv_length') %>% mutate(panel='e', metric='genome-mean log10 inversion length')
)
readr::write_csv(fig2_effects, file.path(OUT_DIR, 'Fig2_effect_sizes_95CI.csv'))
if (file.exists(PGLS_FILE)) {
  pgls_out <- readr::read_csv(PGLS_FILE, show_col_types=FALSE) %>% filter(response %in% c('log10_inv_count_per_Gb','median_log10_inv_length'))
  readr::write_csv(pgls_out, file.path(OUT_DIR, 'Fig2_PGLS_beta_95CI.csv'))
}

save_pdf(Figure2_main, "Fig2_main.pdf", 185, 180)
save_png(Figure2_main, "Fig2_main.png", 185, 180, dpi = 300)

message("Done: Fig2_main.pdf / Fig2_main.png")
