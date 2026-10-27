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

OUT_DIR <- file.path(PROJ_ROOT, 'outputs', 'supplementary')
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

IN_VERIFY <- file.path(PROJ_ROOT, 'allinv.verify.csv')
PER_GENOME <- file.path(PROJ_ROOT, 'per_genome.csv')

per_genome_raw <- load_per_genome(PER_GENOME) 
cols <- c('Metazoa'='#E1812C','Viridiplantae'='#00A087','Fungi'="#3274A1")

############################################################
## Read data
############################################################

VERIFY_dat <- readr::read_csv(IN_VERIFY,show_col_types = FALSE)
validation_genome  <- per_genome_raw %>%
  left_join(VERIFY_dat, by = 'sample_id') 


############################################################
## Prepare data
############################################################

dat <- validation_genome %>%
  mutate(
    length = as.numeric(length),
    status = factor(status,
      levels = c('PASS', 'FAIL')),
    length_class = cut(length,breaks = c(
        -Inf,
        1e3,
        1e4,
        1e5,
        1e6,
        Inf),
      labels = c(
        '≤1 kb',
        '1–10 kb',
        '10–100 kb',
        '100 kb–1 Mb',
        '>1 Mb'),
      right = TRUE))

############################################################
## Summary tables
############################################################

overall <- dat %>%
  group_by(kingdom) %>%
  filter(!is.na(status)) %>%
  count(status, name = 'n') %>%
  mutate(fraction = n / sum(n))

by_length <- dat %>%
  filter( !is.na(status),
    !is.na(length_class)) %>%
  group_by(length_class) %>%
  summarise(
    n_total = n(),
    n_pass = sum(status == 'PASS'),
    n_fail = sum(status == 'FAIL'),
    pass_rate = n_pass / n_total,
    .groups = 'drop')

by_genome <- dat %>%
  filter(!is.na(status)) %>%
  group_by(sample_id,kingdom) %>%
  summarise(n_total = n(),
    n_pass = sum(status == 'PASS'),
    n_fail = sum(status == 'FAIL'),
    pass_rate = n_pass / n_total,
    .groups = 'drop')

############################################################
## A. Overall retention
############################################################

A <- ggplot(overall,
  aes(x = status,
    y = fraction,
    fill = kingdom)) +
  geom_col(width = 0.7,
    position = position_dodge(width = 0.8) ) +
  geom_text(aes(label = paste0(
        sprintf('%.1f%%', 100 * fraction))),
		position = position_dodge(width = 0.8),
    vjust = -0.5,
    size = 2.0) +
  scale_y_continuous(
    limits = c(0, 1),
    labels = scales::percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0.04))) +
  labs(
    x = NULL,
    y = 'Fraction of inversions',
    fill = 'Kingdom') +
	scale_fill_manual(values = cols)+
  theme_nature() +
  theme(legend.position = c(0.75, 0.8),
    legend.direction = "vertical",
    legend.title = element_blank(),
    legend.text = element_text(size = 6) ,
    legend.key.size = grid::unit(3, "mm"))

############################################################
## B. Retention by inversion length
############################################################

B <- ggplot(by_length,
  aes(x = length_class,
    y = pass_rate)) +
  geom_point(size = 2) +
  geom_line(aes(group = 1),
    linewidth = 0.5) +
  scale_y_continuous(
    limits = c(0, 1),
    labels = scales::percent_format(accuracy = 1)) +
  labs( x = NULL,y = 'HiFi validation rate') +
  theme_nature() +
  theme(axis.text.x = element_text(  
    angle =25,                
    hjust = 1,size=5.0),
    legend.position = c(0.8, 0.2),
    legend.direction = "vertical",
    legend.title = element_blank(),
    legend.text = element_text(size = 6) ,
    legend.key.size = grid::unit(3, "mm"))

############################################################
## C. Variation among genomes
############################################################
validation_genome <- dat %>%
  filter(!is.na(status),
    !is.na(kingdom)) %>%
  group_by(sample_id,
    kingdom) %>%
  summarise(n_total = n(),
    n_pass = sum(status == 'PASS'),
    validation_rate = n_pass / n_total,
    .groups = 'drop')%>%
  filter(validation_rate > 0)

C <- ggplot(validation_genome,
  aes(x = kingdom,
    y = validation_rate,
    colour = kingdom)) +
  geom_boxplot(
    width = 0.52,
    outlier.shape = NA,
    colour = 'black') +
  geom_jitter(
    width = 0.10,
    size = 0.8,
    alpha = 0.55) + 
	scale_color_manual(values = cols)+
    scale_y_continuous(limits = c(0, 1),
    labels = scales::percent_format(accuracy = 1)) +
  labs(x = NULL,
    y = 'Per-genome HiFi validation rate') +
  theme_nature() +
  theme(legend.position = 'none')

############################################################
## Assemble
############################################################

SuppFig_validation <-
  (A | B | C) +
  patchwork::plot_annotation(
    tag_levels = 'a',
    theme = theme(
    plot.tag = element_text(
    face = 'bold',
    size = 8)))

save_pdf(SuppFig_validation,
  'SuppFig2_HiFi_validation.pdf',
  180,60)

save_png(SuppFig_validation,
  'SuppFig2_HiFi_validation.png',
  180, 60, dpi = 300)
