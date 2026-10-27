############################################################
## Extended Data Figure 2
## Additional summaries for inversion number and length
##
## Panel order:
##   a  Inversion count per 10,000 genes
##   b  Genome-wide π vs total inversion count
##   c  Proportion of small inversions
##   d  Proportion of medium inversions
##   e  Proportion of large inversions
##   f  Order-level mean inversion count per Gb
##   g  Order-level mean median inversion length
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

per_genome <- load_per_genome()
PGLS_FILE <- file.path(PROJ_ROOT, 'phylogenetic_revision_FINAL', '13_main_Fig2_phylogenetic_correlations.csv')
pgls_b_labels <- read_pgls_labels(PGLS_FILE, 'inv_count_total')
require_columns(
  per_genome,
  c("order", "inv_count_total", "pi_genome_wide", "log10_inv_count_per_Gb", "median_log10_inv_length"),
  "per_genome"
)

if (!all(c("prop_small", "prop_medium", "prop_large") %in% names(per_genome))) {
  stop("Extended Data Fig. 2 requires size-class counts in per_genome.csv to derive prop_small / prop_medium / prop_large.", call. = FALSE)
}

## a
ed2a <- ggplot() + theme_void()
if ("log10_inv_per_10k_genes" %in% names(per_genome)) {
  ed2a_dat <- per_genome %>% filter(!is.na(log10_inv_per_10k_genes))
  ed2a <- boxpanel_kingdom(
    ed2a_dat,
    log10_inv_per_10k_genes,
    ylab = "Inversion count per 10,000 genes",
    x_labels = make_x_labels_with_n(ed2a_dat),
    natural_unit = "count"
  )
} else {
  warning("Column log10_inv_per_10k_genes was not found; panel a is left blank.")
}

## b
ed2b <- corr_panel_kingdom(
  per_genome,
  xvar = "pi_genome_wide",
  yvar = "inv_count_total",
  xlab = "Genome-wide nucleotide diversity (π)",
  ylab = "Total inversion count",
  logy = TRUE,
  method_stats = "pearson",
  legend_pos = "right",
  legend_text_size = 5,
  legend_labels_override = pgls_b_labels,
  natural_y_unit = "count"
)

## c
ed2c_dat <- per_genome %>% filter(!is.na(prop_small))
ed2c <- boxpanel_kingdom(
  ed2c_dat,
  prop_small,
  ylab = "Proportion of small inversions",
  x_labels = make_x_labels_with_n(ed2c_dat)
) +
  scale_y_continuous(labels = percent_format(accuracy = 1))

## d
ed2d_dat <- per_genome %>% filter(!is.na(prop_medium))
ed2d <- boxpanel_kingdom(
  ed2d_dat,
  prop_medium,
  ylab = "Proportion of medium inversions",
  x_labels = make_x_labels_with_n(ed2d_dat)
) +
  scale_y_continuous(labels = percent_format(accuracy = 1))

## e
ed2e_dat <- per_genome %>% filter(!is.na(prop_large))
ed2e <- boxpanel_kingdom(
  ed2e_dat,
  prop_large,
  ylab = "Proportion of large inversions",
  x_labels = make_x_labels_with_n(ed2e_dat)
) +
  scale_y_continuous(labels = percent_format(accuracy = 1))

## order-level summaries
ord_dat <- per_genome %>%
  mutate(order = as.character(order)) %>%
  filter(!is.na(order), order != "") %>%
  group_by(kingdom, order) %>%
  summarise(
    order_mean_log10_inv_perGb = mean(log10_inv_count_per_Gb, na.rm = TRUE),
    order_mean_median_log10_len = mean(median_log10_inv_length, na.rm = TRUE),
    order_mean_prop_large = mean(prop_large, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(kingdom = factor(kingdom, levels = KINGDOM_LEVELS_ALL))

## f
ed2f_dat <- ord_dat %>% filter(!is.na(order_mean_log10_inv_perGb), is.finite(order_mean_log10_inv_perGb))
ed2f <- boxpanel_kingdom(
  ed2f_dat,
  order_mean_log10_inv_perGb,
  ylab = "Inversion count per Gb",
  x_labels = make_x_labels_with_n(ed2f_dat),
  natural_unit = "count"
)

## g
ed2g_dat <- ord_dat %>% filter(!is.na(order_mean_median_log10_len), is.finite(order_mean_median_log10_len))
ed2g <- boxpanel_kingdom(
  ed2g_dat,
  order_mean_median_log10_len,
  ylab = "Median inversion length",
  x_labels = make_x_labels_with_n(ed2g_dat),
  natural_unit = "bp"
)

## h
ed2h_dat <- ord_dat %>% filter(!is.na(order_mean_prop_large), is.finite(order_mean_prop_large))
ed2h_xlab <- make_x_labels_with_n(ed2h_dat)
ed2h <- boxpanel_kingdom(
  ed2h_dat,
  order_mean_prop_large,
  ylab = expression(paste("Proportion of large inversions")),
  x_labels = ed2h_xlab
)

####
row1 <-(ed2a| ed2b|plot_spacer()) + plot_layout(widths = c(1.2, 1.3, 1e-6), guides = "keep")
Extended_Data_Fig2 <-
  row1 /
  (ed2c | ed2d | ed2e) /
  (ed2f | ed2g | ed2h) +
  plot_annotation(tag_levels = "a")

save_pdf(Extended_Data_Fig2, "Extended_Data_Fig2.pdf", 180, 170)
save_png(Extended_Data_Fig2, "Extended_Data_Fig2.png", 180, 170, dpi = 300)

message("Done: Extended_Data_Fig2.pdf / Extended_Data_Fig2.png")
