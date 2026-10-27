get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_flag <- '--file='
  path <- sub(file_flag, '', args[grepl(file_flag, args)])
  if (length(path) > 0) dirname(normalizePath(path[1])) else getwd()
}
script_dir <- get_script_dir()
PROJ_ROOT <- if (basename(script_dir) %in% c('extended_data','supplementary','main_figures')) dirname(script_dir) else script_dir
source(file.path(PROJ_ROOT, 'helpers', 'figure_plot_helpers.R'))
source(file.path(PROJ_ROOT, 'helpers', 'supp_quality_helpers.R'))
source(file.path(PROJ_ROOT, 'helpers', 'figure34_helpers.R'))
OUT_DIR <- file.path(PROJ_ROOT, 'outputs', 'supplementary')
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
IN_PER_GENOME <- file.path(PROJ_ROOT, 'per_genome.csv')
IN_PER_INVERSION <- file.path(PROJ_ROOT, 'per_inversion.csv')
IN_MATING <- file.path(PROJ_ROOT, 'plant_mating.txt')

############################################################
## Mating-system sensitivity
## Remove explicitly selfing plants
############################################################
AGE_THRESHOLD <- 1.0
MAX_POINTS_PER_KINGDOM <- 10000L
SEED <- 123

per_genome_raw <- load_per_genome(IN_PER_GENOME) 
per_inv_raw <- load_per_inversion(IN_PER_INVERSION)
mating <- readr::read_tsv(IN_MATING, show_col_types = FALSE, trim_ws = TRUE) %>%
  transmute(Species = stringr::str_squish(Species),
  mating_system = stringr::str_squish(`Mating system`)) %>%
  distinct(Species, .keep_all = TRUE)

###per_genome dataset
per_genome_mating <- per_genome_raw %>%
  left_join(mating, by = 'Species')

per_genome_outcross  <- per_genome_mating %>%
  filter(!is.na(sample_id),
    kingdom == 'Metazoa' | kingdom == 'Fungi' |
      (kingdom == 'Viridiplantae' &
       mating_system == 'Outcrossing') )

##per_inv_dataset
keep_ids <- per_genome_outcross %>%
  filter( !is.na(sample_id),
          !is.na(pi_genome_wide)) %>%
  pull(sample_id) %>%
  unique()

per_inv_outcross <- per_inv_raw %>%
  filter( sample_id %in% keep_ids,
          kingdom %in% c('Metazoa','Viridiplantae'),
          !is.na(inv_len_bp)) %>%
  mutate( kingdom = factor(kingdom, levels = c('Metazoa','Viridiplantae')),
          inv_len_bp = as.numeric(inv_len_bp),
          log10_len = safe_log10(inv_len_bp),
          pi_ratio = as.numeric(pi_ratio_inv_genome),
          pi_ratio = if_else( !is.finite(pi_ratio) | pi_ratio <= 0,
                              NA_real_,
                              pi_ratio),
          age_class = case_when(is.na(pi_ratio) ~ NA_character_,
                                pi_ratio <= AGE_THRESHOLD ~ 'Young',
                                pi_ratio > AGE_THRESHOLD ~ 'Old'),
          age_class = factor(age_class, levels = c('Young','Old')),
          size_class = case_when( inv_len_bp < 1e4 ~ 'Small',
                                  inv_len_bp >= 1e4 & inv_len_bp <= 1e6 ~ 'Medium',
                                  inv_len_bp > 1e6 ~ 'Large',
                                  TRUE ~ NA_character_),
          size_class = factor(size_class, levels = c('Small','Medium','Large'))  )
## explicitly annotated as Selfing. Unmatched/other plants are retained.
per_genome_noself <- per_genome_mating %>%
  filter(!(kingdom == 'Viridiplantae' & mating_system == 'Selfing'))

noself_ids <- per_genome_noself$sample_id
per_inv_noself <- per_inv_raw %>%
  filter(sample_id %in% noself_ids,
         kingdom %in% c('Metazoa','Viridiplantae'),
         !is.na(inv_len_bp)) %>%
  mutate(kingdom = factor(kingdom, levels = c('Metazoa','Viridiplantae')),
         inv_len_bp = as.numeric(inv_len_bp),
         log10_len = safe_log10(inv_len_bp),
         pi_ratio = as.numeric(pi_ratio_inv_genome),
         pi_ratio = if_else(!is.finite(pi_ratio) | pi_ratio <= 0, NA_real_, pi_ratio),
         age_class = case_when(is.na(pi_ratio) ~ NA_character_,
                               pi_ratio <= AGE_THRESHOLD ~ 'Low',
                               pi_ratio > AGE_THRESHOLD ~ 'Elevated'),
         age_class = factor(age_class, levels = c('Low','Elevated')),
         size_class = case_when(inv_len_bp < 1e4 ~ 'Small',
                                inv_len_bp >= 1e4 & inv_len_bp <= 1e6 ~ 'Medium',
                                inv_len_bp > 1e6 ~ 'Large',
                                TRUE ~ NA_character_),
         size_class = factor(size_class, levels = c('Small','Medium','Large')))
############################################################
## Mating-system diagnostic
############################################################
diag_dat <- per_genome_mating %>%
  filter( kingdom == 'Viridiplantae',
  mating_system %in% c('Outcrossing','Mixed','Selfing')) %>%
  mutate(mating_system = factor(mating_system,
         levels = c('Outcrossing','Mixed','Selfing')))

###
mating_cols <- c("Outcrossing" = "#0B6623",
  "Mixed"       = "#7BB661",
  "Selfing"     = "#E6A15C")

###plot function
boxpanel_mating <- function(df, y, ylab, x_labels, add_p = TRUE, natural_unit = NULL) {
  df2 <- df %>% filter(!is.na({{ y }}), is.finite({{ y }}))
  ycol <- rlang::as_name(rlang::enquo(y))
  
  ptab <- manual_pairwise_wilcox(df2, ycol) %>%
    add_y_positions(max(df2[[ycol]], na.rm = TRUE))
  
  p <- ggplot(df2, aes(x = mating_system, y = {{ y }}, fill = mating_system, color = mating_system)) +
    geom_boxplot(width = 0.6, alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.15, size = 1.2, alpha = 0.7) +
    scale_fill_manual(values = mating_cols) +
    scale_color_manual(values = mating_cols) +
    scale_x_discrete(labels = x_labels) +
    labs(x = NULL, y = ylab) +
    theme_nature()
  
  if (add_p) p <- maybe_add_pvalues(p, ptab)
  if (!is.null(natural_unit)) p <- backtransform_log_axis(p, df2[[ycol]], axis="y", unit=natural_unit, max_labels=4)
  p
}

###polting
DiagA <- boxpanel_mating(diag_dat,pi_genome_wide,
  ylab = expression(pi[genome]),
  x_labels = make_x_labels_with_n(diag_dat))


DiagB <- boxpanel_mating(diag_dat,log10_inv_count_per_Gb,
  ylab = "Inversion count per Gb",
  x_labels = make_x_labels_with_n(diag_dat), natural_unit="count")

DiagC <- boxpanel_mating(diag_dat,median_log10_inv_length,
  ylab = "Median inversion length",
  x_labels = make_x_labels_with_n(diag_dat), natural_unit="bp")

MatingDiag <- (DiagA | DiagB | DiagC) +
  plot_annotation(tag_levels = 'a',
    theme = theme(plot.tag = element_text(face = 'bold', size = 8)))

save_pdf(MatingDiag, 'SuppFig12_MatingSystem_diagnostic.pdf', 180, 60)
save_png(MatingDiag, 'SuppFig12_MatingSystem_diagnostic.png', 180, 60, dpi = 300)

diag_summary <- diag_dat %>%
  group_by(mating_system) %>%
  summarise(n = n(),
    median_pi_genome = median(pi_genome_wide, na.rm = TRUE),
    median_inv_count_per_Gb = median(inv_count_per_Gb, na.rm = TRUE),
    median_inv_len_bp = median(inv_len_median_bp, na.rm = TRUE),
    .groups = 'drop')

readr::write_csv(diag_summary,
  file.path(OUT_DIR, 'MatingSystem_diagnostic_summary.csv'))

############################################################
## Mating-system repeat of Supplementary Figure 2
## Complete Figure 2-style analysis, identical 3 x 3 layout.
############################################################
ds_ms_noself <- build_quality_dataset(per_genome_noself,
  per_inv_noself)
build_quality_figure( ds_ms_noself,'SuppFig_Mating_noSelfing_Fig2')
