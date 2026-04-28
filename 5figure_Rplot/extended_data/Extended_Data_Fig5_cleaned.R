
get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_flag <- '--file='
  path <- sub(file_flag, '', args[grepl(file_flag, args)])
  if (length(path) > 0) dirname(normalizePath(path[1])) else getwd()
}
SCRIPT_DIR <- get_script_dir()
PROJ_ROOT <- normalizePath(file.path(SCRIPT_DIR, '..'), winslash = '/', mustWork = FALSE)
source(file.path(PROJ_ROOT, 'helpers', 'figure_plot_helpers.R'))
source(file.path(PROJ_ROOT, 'helpers', 'figure34_helpers.R'))
OUT_DIR <- file.path(PROJ_ROOT, 'outputs', 'extended_data')
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
IN_PER_GENOME <- file.path(PROJ_ROOT, 'per_genome.csv')
IN_PER_INVERSION <- file.path(PROJ_ROOT, 'per_inversion.csv')


############################################################
## Extended Data Figure 5
############################################################
IN_KAKS_SHIFT <- file.path(PROJ_ROOT, 'inversion_kaks_shift.csv')
if (!file.exists(IN_KAKS_SHIFT)) stop('Extended Data Fig. 5 expects inversion_kaks_shift.csv in the package root.', call. = FALSE)
per_inv <- load_per_inversion_age(IN_PER_INVERSION)
len_gene <- summarise_by_lenbin(per_inv %>% filter(!is.na(log10_genes)), log10_genes)
ed5a <- plot_lenbin_trait(len_gene, ylab = expression('Number of genes per inversion (log'[10] * ')'))
ed5b <- plot_lenclass_box(per_inv, te_prop_union, ylab = 'TE proportion')
inv_age_scaling_w <- per_inv %>% filter(!is.na(pi_ratio_inv_flank), pi_ratio_inv_flank > 0) %>% winsorize_by_group_bin(pi_ratio_inv_flank, group_col = 'kingdom', bin_col = pi_ratio, p = 0.999)
age_pi <- summarise_by_agebin(inv_age_scaling_w, pi_ratio_inv_flank)
ed5c <- plot_agebin_trait(age_pi, ylab = expression(pi[inversion] / pi[flank]), hline1 = 1)
kaks_shift <- readr::read_tsv(IN_KAKS_SHIFT, show_col_types = FALSE) %>% mutate(kaks_shift_ratio_median = suppressWarnings(as.numeric(kaks_shift_ratio_median)), n_genes_found = suppressWarnings(as.numeric(n_genes_found)))
inv_meta_shift <- per_inv %>% transmute(sample_id, inv_id, kingdom, inv_len_bp, size_bin = case_when(inv_len_bp < 1e5 ~ 'small', inv_len_bp < 1e6 ~ 'medium', TRUE ~ 'large'), size_bin = factor(size_bin, levels = c('small','medium','large')))
shift_df <- kaks_shift %>% inner_join(inv_meta_shift, by = c('sample_id','inv_id')) %>% filter(kingdom %in% c('Metazoa','Viridiplantae'), status == 'ok', !is.na(size_bin), !is.na(n_genes_found), n_genes_found >= 1)
ed5d <- plot_density_kingdom_ratio(shift_df, 'kaks_shift_ratio_median', xlab = expression(paste('Relative ', omega, ' of inversion genes')), xlim = c(0,2))
ed5e <- plot_density_size(shift_df, 'kaks_shift_ratio_median', xlab = expression(paste('Relative ', omega, ' of inversion genes')), xlim = c(0,2), facet_by_kingdom = TRUE)
row2 <- ed5d + ed5e + plot_layout(widths = c(1,2))
ExtendedDataFig5 <- (ed5a | ed5b | ed5c) / row2 + plot_annotation(tag_levels = 'a', theme = theme(plot.tag = element_text(face = 'bold', size = 8)))
save_pdf(ExtendedDataFig5, 'Extended_Data_Fig5.pdf', 180, 120); save_png(ExtendedDataFig5, 'Extended_Data_Fig5.png', 180, 120, dpi = 300); writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, 'Extended_Data_Fig5_sessionInfo.txt'))
