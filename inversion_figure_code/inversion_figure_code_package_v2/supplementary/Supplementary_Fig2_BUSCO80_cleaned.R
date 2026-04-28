
get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_flag <- '--file='
  path <- sub(file_flag, '', args[grepl(file_flag, args)])
  if (length(path) > 0) dirname(normalizePath(path[1])) else getwd()
}
SCRIPT_DIR <- get_script_dir()
PROJ_ROOT <- normalizePath(file.path(SCRIPT_DIR, '..'), winslash = '/', mustWork = FALSE)
source(file.path(PROJ_ROOT, 'helpers', 'figure_plot_helpers.R'))
source(file.path(PROJ_ROOT, 'helpers', 'supp_quality_helpers.R'))
OUT_DIR <- file.path(PROJ_ROOT, 'outputs', 'supplementary')
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
IN_PER_GENOME <- file.path(PROJ_ROOT, 'per_genome.csv')
IN_PER_INVERSION <- file.path(PROJ_ROOT, 'per_inversion.csv')


############################################################
## SuppFig2_BUSCO80
############################################################
per_genome_raw <- load_per_genome(IN_PER_GENOME) %>% mutate(BUSCO_completeness = readr::parse_number(as.character(BUSCO_completeness)), BUSCO_duplication = readr::parse_number(as.character(BUSCO_duplication)))
per_inv_raw <- load_per_inversion(IN_PER_INVERSION)
ds <- build_quality_dataset(per_genome_raw, per_inv_raw, mode = 'BUSCO80')
readr::write_csv(ds$per_genome, file.path(OUT_DIR, 'SuppFig2_BUSCO80_per_genome.csv'))
readr::write_csv(ds$per_inv, file.path(OUT_DIR, 'SuppFig2_BUSCO80_per_inversion.csv'))
build_quality_figure(ds, 'SuppFig2_BUSCO80')
writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, 'SuppFig2_BUSCO80_sessionInfo.txt'))
