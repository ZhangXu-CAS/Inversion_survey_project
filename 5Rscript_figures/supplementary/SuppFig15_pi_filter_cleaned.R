get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_flag <- '--file='
  path <- sub(file_flag, '', args[grepl(file_flag, args)])
  if (length(path) > 0) dirname(normalizePath(path[1])) else getwd()
}

script_dir <- get_script_dir()
PROJ_ROOT <- if (basename(script_dir) %in% c('extended_data','supplementary','main_figures')) dirname(script_dir) else script_dir

source(file.path(PROJ_ROOT, 'helpers', 'figure_plot_helpers.R'))
source(file.path(PROJ_ROOT, 'helpers', 'figure34_helpers.R'))

OUT_DIR <- file.path(PROJ_ROOT, 'outputs', 'supplementary')
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

IN_PER_GENOME <- file.path(PROJ_ROOT, 'per_genome.csv')
IN_PER_INVERSION <- file.path(PROJ_ROOT, 'per_inversion.csv')
IN_MATING <- file.path(PROJ_ROOT, 'plant_mating.txt')


############################################################
## Supplementary Figure 5
## Robustness to genome-wide haplotype-divergence thresholds
## and plant mating system
############################################################

PI_GENOME_MIN4 <- 1e-4
PI_GENOME_MIN3 <- 1e-3


############################################################
## Read data
############################################################

per_inv_raw <- readr::read_csv(
  IN_PER_INVERSION,
  show_col_types = FALSE
)

per_genome_raw <- readr::read_csv(
  IN_PER_GENOME,
  show_col_types = FALSE
)

mating <- readr::read_tsv(
  IN_MATING,
  show_col_types = FALSE,
  trim_ws = TRUE
) %>%
  transmute(
    Species = stringr::str_squish(Species),
    mating_system = stringr::str_squish(`Mating system`)
  ) %>%
  distinct(Species, .keep_all = TRUE)


############################################################
##pi-filtered datasets
############################################################

per_genome_filt4 <- per_genome_raw %>%
  mutate(
    pi_genome_wide = as.numeric(pi_genome_wide)
  ) %>%
  filter(
    !is.na(sample_id),
    !is.na(pi_genome_wide),
    pi_genome_wide >= PI_GENOME_MIN4
  )

per_genome_filt3 <- per_genome_raw %>%
  mutate(
    pi_genome_wide = as.numeric(pi_genome_wide)
  ) %>%
  filter(
    !is.na(sample_id),
    !is.na(pi_genome_wide),
    pi_genome_wide >= PI_GENOME_MIN3
  )


############################################################
## Strict-outcrossing mating-system dataset
############################################################

per_genome_mating <- per_genome_raw %>%
  mutate(
    Species = stringr::str_squish(Species),
    pi_genome_wide = as.numeric(pi_genome_wide)
  ) %>%
  left_join(
    mating,
    by = 'Species'
  )

per_genome_outcross <- per_genome_mating %>%
  filter(
    !is.na(sample_id),
    !is.na(pi_genome_wide),
    kingdom == 'Metazoa' | kingdom == 'Fungi' |
      (  kingdom == 'Viridiplantae' &
          mating_system == 'Outcrossing') )

############################################################
## Row 1
## All taxa with pi >= 1e-4
############################################################

S0 <- corr_panel_kingdom(
  per_genome_filt4,
  xvar = 'pi_genome_wide',
  yvar = 'inv_count_per_Gb',
  xlab = 'Genome-wide haplotype divergence (π)',
  ylab = 'Inversion count per Gb',
  logy = TRUE,
  method_stats = 'pearson',
  legend_pos = c(0.70, 0.20),
  legend_text_size = 5,
  natural_y_unit = 'count'
)

S1 <- corr_panel_kingdom(
  per_genome_filt4,
  xvar = 'pi_genome_wide',
  yvar = 'inv_count_total',
  xlab = 'Genome-wide haplotype divergence (π)',
  ylab = 'Total inversion count',
  logy = TRUE,
  method_stats = 'pearson',
  legend_pos = c(0.70, 0.20),
  legend_text_size = 5,
  natural_y_unit = 'count'
)

S2 <- corr_panel_kingdom(
  per_genome_filt4,
  xvar = 'pi_genome_wide',
  yvar = 'median_log10_inv_length',
  xlab = 'Genome-wide haplotype divergence (π)',
  ylab = 'Median inversion length',
  logy = FALSE,
  legend_pos = c(0.70, 0.85),
  legend_text_size = 5,
  natural_y_unit = 'bp'
)


############################################################
## Row 2
## All taxa with pi >= 1e-3
############################################################

S3 <- corr_panel_kingdom(
  per_genome_filt3,
  xvar = 'pi_genome_wide',
  yvar = 'inv_count_per_Gb',
  xlab = 'Genome-wide haplotype divergence (π)',
  ylab = 'Inversion count per Gb',
  logy = TRUE,
  method_stats = 'pearson',
  legend_pos = c(0.70, 0.20),
  legend_text_size = 5,
  natural_y_unit = 'count'
)

S4 <- corr_panel_kingdom(
  per_genome_filt3,
  xvar = 'pi_genome_wide',
  yvar = 'inv_count_total',
  xlab = 'Genome-wide haplotype divergence (π)',
  ylab = 'Total inversion count',
  logy = TRUE,
  method_stats = 'pearson',
  legend_pos = c(0.70, 0.20),
  legend_text_size = 5,
  natural_y_unit = 'count'
)

S5 <- corr_panel_kingdom(
  per_genome_filt3,
  xvar = 'pi_genome_wide',
  yvar = 'median_log10_inv_length',
  xlab = 'Genome-wide haplotype divergence (π)',
  ylab = 'Median inversion length',
  logy = FALSE,
  legend_pos = c(0.70, 0.85),
  legend_text_size = 5,
  natural_y_unit = 'bp'
)


############################################################
## Row 3
## Strict mating-system control
##
## Plants = Outcrossing only
## Animals = all Metazoa
############################################################

S6 <- corr_panel_kingdom(
  per_genome_outcross,
  xvar = 'pi_genome_wide',
  yvar = 'inv_count_per_Gb',
  xlab = 'Genome-wide haplotype divergence (π)',
  ylab = 'Inversion count per Gb',
  logy = TRUE,
  method_stats = 'pearson',
  legend_pos = c(0.70, 0.20),
  legend_text_size = 5,
  natural_y_unit = 'count'
)

S7 <- corr_panel_kingdom(
  per_genome_outcross,
  xvar = 'pi_genome_wide',
  yvar = 'inv_count_total',
  xlab = 'Genome-wide haplotype divergence (π)',
  ylab = 'Total inversion count',
  logy = TRUE,
  method_stats = 'pearson',
  legend_pos = c(0.70, 0.20),
  legend_text_size = 5,
  natural_y_unit = 'count'
)

S8 <- corr_panel_kingdom(
  per_genome_outcross,
  xvar = 'pi_genome_wide',
  yvar = 'median_log10_inv_length',
  xlab = 'Genome-wide haplotype divergence (π)',
  ylab = 'Median inversion length',
  logy = FALSE,
  legend_pos = c(0.70, 0.85),
  legend_text_size = 5,
  natural_y_unit = 'bp'
)


############################################################
## Assemble Supplementary Figure 5
############################################################

SuppFig15_pi_filter <-
  (S0 | S1 | S2) /
  (S3 | S4 | S5) /
  (S6 | S7 | S8) +
  plot_annotation(
    tag_levels = 'a',
    theme = theme(
      plot.tag = element_text(
        face = 'bold',
        size = 8)))

save_pdf( SuppFig14_pi_filter,
  'SuppFig15_pi_filter.pdf',180,180)

save_png( SuppFig15_pi_filter,'SuppFig15_pi_filter.png',180,180,dpi = 300)

writeLines( capture.output(sessionInfo()),
  file.path(OUT_DIR,'SuppFig15_pi_filter_sessionInfo.txt'))