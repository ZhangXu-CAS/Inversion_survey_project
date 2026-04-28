
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
OUT_DIR <- file.path(PROJ_ROOT, 'outputs', 'supplementary')
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
IN_PER_GENOME <- file.path(PROJ_ROOT, 'per_genome.csv')
IN_PER_INVERSION <- file.path(PROJ_ROOT, 'per_inversion.csv')


############################################################
## Supplementary Figure 5
############################################################
PI_GENOME_MIN4 <- 1e-4; PI_GENOME_MIN3 <- 1e-3
per_inv_raw <- readr::read_csv(IN_PER_INVERSION, show_col_types = FALSE); per_genome_raw <- readr::read_csv(IN_PER_GENOME, show_col_types = FALSE)
id_candidates <- c('sample_id','genome_id','genome','assembly_id','assembly','accession','acc','id'); id_col <- intersect(id_candidates, intersect(names(per_genome_raw), names(per_inv_raw)))[1]; if (is.na(id_col)) stop('Cannot find a shared ID column between per_genome.csv and per_inversion.csv.', call. = FALSE)
pi_candidates <- c('pi_genome_wide','pi_genomewide','pi_genome','pi_genome_all','pi_genomewide_mean'); pi_col <- intersect(pi_candidates, names(per_genome_raw))[1]; if (is.na(pi_col)) stop('Cannot find a genome-wide pi column in per_genome.csv.', call. = FALSE)
per_genome_filt4 <- per_genome_raw %>% mutate(pi_genome_wide = as.numeric(.data[[pi_col]])) %>% filter(!is.na(.data[[id_col]]), !is.na(pi_genome_wide), pi_genome_wide >= PI_GENOME_MIN4)
per_genome_filt3 <- per_genome_raw %>% mutate(pi_genome_wide = as.numeric(.data[[pi_col]])) %>% filter(!is.na(.data[[id_col]]), !is.na(pi_genome_wide), pi_genome_wide >= PI_GENOME_MIN3)
readr::write_csv(per_genome_filt4, file.path(OUT_DIR, 'SuppFig5_per_genome_pi_ge_1e4.csv')); readr::write_csv(per_genome_filt3, file.path(OUT_DIR, 'SuppFig5_per_genome_pi_ge_1e3.csv'))
S0 <- corr_panel_kingdom(per_genome_filt4, xvar = 'pi_genome_wide', yvar = 'inv_count_per_Gb', xlab = 'Genome-wide nucleotide diversity (π)', ylab = expression(paste('Inversion count per Gb (log'[10], ')')), logy = TRUE, method_stats = 'pearson', legend_pos = c(0.70,0.20), legend_text_size = 5)
S1 <- corr_panel_kingdom(per_genome_filt4, xvar = 'pi_genome_wide', yvar = 'inv_count_total', xlab = 'Genome-wide nucleotide diversity (π)', ylab = expression(paste('Total inversion count (log'[10], ')')), logy = TRUE, legend_pos = c(0.70,0.20), legend_text_size = 5)
S2 <- corr_panel_kingdom(per_genome_filt4, xvar = 'pi_genome_wide', yvar = 'median_log10_inv_length', xlab = 'Genome-wide nucleotide diversity (π)', ylab = expression(paste('Median inversion length (log'[10], ' bp)')), logy = FALSE, legend_pos = c(0.70,0.85), legend_text_size = 5)
S3 <- corr_panel_kingdom(per_genome_filt3, xvar = 'pi_genome_wide', yvar = 'inv_count_per_Gb', xlab = 'Genome-wide nucleotide diversity (π)', ylab = expression(paste('Inversion count per Gb (log'[10], ')')), logy = TRUE, method_stats = 'pearson', legend_pos = c(0.70,0.20), legend_text_size = 5)
S4 <- corr_panel_kingdom(per_genome_filt3, xvar = 'pi_genome_wide', yvar = 'inv_count_total', xlab = 'Genome-wide nucleotide diversity (π)', ylab = expression(paste('Total inversion count (log'[10], ')')), logy = TRUE, legend_pos = c(0.70,0.20), legend_text_size = 5)
S5 <- corr_panel_kingdom(per_genome_filt3, xvar = 'pi_genome_wide', yvar = 'median_log10_inv_length', xlab = 'Genome-wide nucleotide diversity (π)', ylab = expression(paste('Median inversion length (log'[10], ' bp)')), logy = FALSE, legend_pos = c(0.70,0.85), legend_text_size = 5)
SuppFig5 <- (S0 | S1 | S2) / (S3 | S4 | S5) + plot_annotation(tag_levels = 'a', theme = theme(plot.tag = element_text(face = 'bold', size = 8)))
save_pdf(SuppFig5, 'SuppFig5.pdf', 180, 120); save_png(SuppFig5, 'SuppFig5.png', 180, 120, dpi = 300); writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, 'SuppFig5_sessionInfo.txt'))
