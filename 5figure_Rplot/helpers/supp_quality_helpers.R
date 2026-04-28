
suppressPackageStartupMessages({
  library(tidyverse)
  library(ggpubr)
  library(patchwork)
  library(scales)
})

make_x_labels_n <- function(df) {
  df %>% count(kingdom, name = 'n') %>% mutate(label = paste0(as.character(kingdom), '\n(', n, ')')) %>% select(kingdom, label) %>% deframe()
}

panel_box_quality <- function(df, y, ylab) {
  df2 <- df %>% filter(!is.na(.data[[y]]), is.finite(.data[[y]]))
  if (nrow(df2) == 0) return(ggplot() + theme_void())
  xlab <- make_x_labels_n(df2)
  ptab <- manual_pairwise_wilcox(df2, y) %>% add_y_positions(max(df2[[y]], na.rm = TRUE))
  p <- ggplot(df2, aes(x = kingdom, y = .data[[y]], fill = kingdom, color = kingdom)) +
    geom_boxplot(width = 0.6, alpha = 0.7, outlier.shape = NA) + geom_jitter(width = 0.15, size = 1.2, alpha = 0.7) +
    scale_fill_manual(values = kingdom_cols) + scale_color_manual(values = kingdom_cols) + scale_x_discrete(labels = xlab) + labs(x = NULL, y = ylab) + theme_nature()
  maybe_add_pvalues(p, ptab)
}

panel_prop_quality <- function(df, y, ylab) panel_box_quality(df, y, ylab) + scale_y_continuous(labels = percent_format(accuracy = 1))

panel_density_len_quality <- function(per_inv, add_legend = TRUE) {
  p <- ggplot(per_inv, aes(x = log10_len, color = kingdom, fill = kingdom)) + geom_density(alpha = 0.25, adjust = 1.0) + scale_color_manual(values = kingdom_cols) + scale_fill_manual(values = kingdom_cols) + labs(x = expression(paste('Inversion length (log'[10], ' bp)')), y = 'Density') + theme_nature()
  if (add_legend) p <- p + theme(legend.position = c(0.98,0.98), legend.justification = c(1,1), legend.direction='vertical', legend.title = element_blank(), legend.text = element_text(size=6), legend.key.size = grid::unit(3,'mm'))
  p
}

build_quality_dataset <- function(per_genome_raw, per_inv_raw, mode = c('BUSCO80','BUSCO80_DUP50','LOCO_noInsects')) {
  mode <- match.arg(mode); pg <- per_genome_raw
  if (mode %in% c('BUSCO80','BUSCO80_DUP50')) pg <- pg %>% filter(!is.na(BUSCO_completeness), BUSCO_completeness >= 80)
  if (mode == 'BUSCO80_DUP50') pg <- pg %>% filter(!is.na(BUSCO_duplication), BUSCO_duplication <= 50)
  if (mode == 'LOCO_noInsects') {
    set.seed(20260127)
    insects_ids <- pg %>% filter(kingdom == 'Metazoa', !is.na(major_clade), major_clade == 'Insects') %>% distinct(sample_id) %>% pull(sample_id)
    n_keep <- min(50L, length(insects_ids)); keep <- if (n_keep > 0) sample(insects_ids, size = n_keep, replace = FALSE) else character(0)
    pg <- pg %>% filter(!(kingdom == 'Metazoa' & !is.na(major_clade) & major_clade == 'Insects') | sample_id %in% keep)
  }
  pg <- derive_size_props(pg) %>% mutate(log10_inv_count_per_Gb = if_else(is.finite(inv_count_per_Gb) & inv_count_per_Gb > 0, safe_log10(inv_count_per_Gb), NA_real_), log10_inv_count_total = if_else(is.finite(inv_count_total) & inv_count_total > 0, safe_log10(inv_count_total), NA_real_))
  pi <- per_inv_raw %>% filter(sample_id %in% unique(pg$sample_id))
  list(per_genome = pg, per_inv = pi)
}

build_quality_figure <- function(ds, filename_stub) {
  pg <- ds$per_genome; pi <- ds$per_inv
  p1 <- panel_box_quality(pg, 'log10_inv_count_per_Gb', expression(paste('Inversion count per Gb (log'[10], ')')))
  p2 <- panel_box_quality(pg, 'log10_inv_count_total', expression(paste('Total inversion count (log'[10], ')')))
  p3 <- if ('log10_inv_per_10k_genes' %in% names(pg)) panel_box_quality(pg, 'log10_inv_per_10k_genes', expression(paste('Inversion count per 10,000 genes (log'[10], ')'))) else ggplot() + theme_void()
  p4 <- panel_box_quality(pg, 'median_log10_inv_length', expression(paste('Median inversion length (log'[10], ' bp)')))
  p5 <- panel_box_quality(pg, 'mean_log10_inv_length', expression(paste('Mean inversion length (log'[10], ' bp)')))
  p6 <- panel_density_len_quality(pi, add_legend = TRUE)
  p7 <- panel_prop_quality(pg, 'prop_small', 'Proportion of small inversions')
  p8 <- panel_prop_quality(pg, 'prop_medium', 'Proportion of medium inversions')
  p9 <- panel_prop_quality(pg, 'prop_large', 'Proportion of large inversions')
  fig <- (p1 | p2 | p3) / (p4 | p5 | p6) / (p7 | p8 | p9) + plot_annotation(tag_levels = 'a', theme = theme(plot.tag = element_text(face = 'bold', size = 8)))
  save_pdf(fig, paste0(filename_stub, '.pdf'), 180, 170); save_png(fig, paste0(filename_stub, '.png'), 180, 170, dpi = 300); invisible(fig)
}
