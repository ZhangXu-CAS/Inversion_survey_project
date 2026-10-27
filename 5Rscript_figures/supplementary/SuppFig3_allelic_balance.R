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

IN_INV <- file.path(
  PROJ_ROOT,
  'allelic_balance_inversions.tsv'
)

IN_SAMPLE <- file.path(
  PROJ_ROOT,
  'allelic_balance_samples.tsv'
)

cols <- c('Metazoa'='#E1812C','Viridiplantae'='#00A087','Fungi'="#3274A1")
############################################################
## Read
############################################################

inv <- readr::read_tsv(
  IN_INV,
  show_col_types = FALSE
)

samples <- readr::read_tsv(
  IN_SAMPLE,
  show_col_types = FALSE
)


############################################################
## Required columns
############################################################

required_inv <- c(
  'sample_id',
  'kingdom',
  'length_bp',
  'length_class',
  'AB_MAPQ20',
  'centered_AB_MAPQ20',
  'AB_evaluable_10_MAPQ20',
  'H1_mean_depth',
  'H2_mean_depth',
  'depth_ratio'
)

missing_inv <- setdiff(
  required_inv,
  names(inv)
)

if (length(missing_inv) > 0) {
  stop(
    paste(
      'Missing inversion columns:',
      paste(missing_inv, collapse = ', ')
    ),
    call. = FALSE
  )
}


############################################################
## Primary analysis set
############################################################
length_levels <- c(
  '<=1 kb',
  '1–10 kb',
  '10–100 kb',
  '100 kb–1 Mb',
  '>1 Mb'
)

length_labels <- c(
  '≤1kb',
  '1–10kb',
  '10–100kb',
  '100kb–1Mb',
  '>1Mb'
)


dat <- inv %>%
  filter(
    AB_evaluable_10_MAPQ20,
    !is.na(AB_MAPQ20),AB_MAPQ20 != 1.00,
    !is.na(centered_AB_MAPQ20)
  ) %>%
  mutate(
    length_class = factor(
      length_class,
      levels = length_levels,
      labels = length_labels
    ),
    
    abs_centered_AB =
      abs(centered_AB_MAPQ20),
    
    log2_depth_ratio =
      log2(depth_ratio)
  )
############################################################
## Genome-level summaries
############################################################

genome_summary <- dat %>%
  group_by(
    sample_id,
    kingdom
  ) %>%
  summarise(
    median_AB =
      median(AB_MAPQ20, na.rm = TRUE),

    median_centered_AB =
      median(centered_AB_MAPQ20, na.rm = TRUE),

    median_abs_centered_AB =
      median(
        abs(centered_AB_MAPQ20),
        na.rm = TRUE
      ),
    median_log2_depth_ratio =
      median(
        log2_depth_ratio[
          is.finite(log2_depth_ratio)
        ],
        na.rm = TRUE
      ),
      .groups = 'drop'
  )


############################################################
## A. Genome-wide competitive haplotype representation
############################################################

A <- ggplot(
  samples,
  aes(
    x = kingdom,
    y = genome_AB,
    colour = kingdom
  )
) +
  geom_hline(
    yintercept = 0.5,
    linetype = 2,
    linewidth = 0.5
  ) +
  geom_boxplot(
    width = 0.52,
    outlier.shape = NA,
    colour = 'black'
  ) +scale_color_manual(values = cols)+
  geom_jitter(
    width = 0.10,
    size = 1.3,
    alpha = 0.75
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = c(0, 0.25, 0.5, 0.75, 1)
  ) +
  labs(
    x = NULL,
    y = 'Genome-wide H1 allelic balance'
  ) +
  theme_nature() +
  theme(
    legend.position = 'none'
  )


############################################################
## B. Inversion-level allelic balance
############################################################

B <- ggplot(
  dat,
  aes(x = AB_MAPQ20)
) +
  geom_histogram(
    binwidth = 0.025,
    boundary = 0
  ) +
  geom_vline(
    xintercept = 0.5,
    linetype = 2,
    linewidth = 0.5
  ) +
  scale_x_continuous(
    limits = c(0, 1),
    breaks = c(0, 0.25, 0.5, 0.75, 1)
  ) +
  labs(
    x = 'Inversion H1 allelic balance',
    y = 'Number of inversions'
  ) +
  theme_nature()


############################################################
## C. Inversion deviation from genome-wide expectation
############################################################

C <- ggplot(
  dat,
  aes(x = centered_AB_MAPQ20)
) +
  geom_histogram(
    binwidth = 0.025,
    boundary = 0
  ) +
  geom_vline(
    xintercept = 0,
    linetype = 2,
    linewidth = 0.5
  ) +
  labs(
    x = 'Inversion AB − genome-wide AB',
    y = 'Number of inversions'
  ) +
  theme_nature()


############################################################
## D. Genome-centred deviation by length and kingdom
############################################################

length_summary <- dat %>%
  filter(!is.na(length_class)) %>%
  group_by(
    sample_id,
    length_class
  ) %>%
  summarise(
    median_abs_centered_AB =
      median(
        abs(centered_AB_MAPQ20),
        na.rm = TRUE
      ),
    .groups = 'drop'
  )

D <- ggplot(
  length_summary,
  aes(
    x = length_class,
    y = median_abs_centered_AB
  )
) +
  geom_boxplot(
    aes(
      group = interaction(
        length_class
      )
    ),
    position = position_dodge(
      width = 0.72
    ),
    width = 0.55,
    outlier.shape = NA
  ) +
  geom_point(
    position = position_jitterdodge(
      jitter.width = 0.08,
      dodge.width = 0.72
    ),
    size = 0.65,
    alpha = 0.55
  ) +
  labs(
    x = NULL,
    y = 'Per-genome median\n|inversion AB − genome AB|'
  ) +
  theme_nature() +
  theme(
    axis.text.x = element_text(
      angle = 15,
      hjust = 1))


############################################################
## Assemble
############################################################

SuppFig_AB <-
  (A | B) /
  (C | D) +
  patchwork::plot_annotation(
    tag_levels = 'a',
    theme = theme(
      plot.tag = element_text(
        face = 'bold',
        size = 8
      )
    )
  )

save_pdf(
  SuppFig_AB,
  'SuppFig3_allelic_balance.pdf',
  120,
  110
)

save_png(
  SuppFig_AB,
  'SuppFig3_allelic_balance.png',
  120,
  110,
  dpi = 300
)