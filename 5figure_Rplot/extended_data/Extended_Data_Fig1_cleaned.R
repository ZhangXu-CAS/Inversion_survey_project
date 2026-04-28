############################################################
## Extended Data Figure 1
## Assembly quality and genome architecture summaries
##
## Panel order:
##   a  Assembly size vs scaffold N50
##   b  BUSCO completeness
##   c  TE proportion
##   d  Gene density distribution
##   e  Genome-wide nucleotide diversity (π)
##   f  Genome-wide indel rate
##
## Note:
##   Main Figure 1 is Illustrator-drawn and is not included here.
############################################################


get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_flag <- '--file='
  path <- sub(file_flag, '', args[grepl(file_flag, args)])
  if (length(path) > 0) dirname(normalizePath(path[1])) else getwd()
}
SCRIPT_DIR <- get_script_dir()
PROJ_ROOT <- normalizePath(file.path(SCRIPT_DIR, '..'), winslash = '/', mustWork = FALSE)
source(file.path(PROJ_ROOT, 'helpers', 'figure_plot_helpers.R'))
OUT_DIR <- file.path(PROJ_ROOT, 'outputs', 'extended_data')
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

per_genome <- load_per_genome()
require_columns(
  per_genome,
  c("asm_total_len_bp", "asm_scaffold_N50_Mb", "BUSCO_completeness", "TE_total_pct", "gene_density", "pi_genome_wide", "indel_rate"),
  "per_genome"
)

p_label <- function(p) {
  stars <- dplyr::case_when(
    p < 1e-3 ~ "***",
    p < 1e-2 ~ "**",
    p < 5e-2 ~ "*",
    TRUE ~ ""
  )

  dplyr::if_else(
    p < 1e-3,
    "p < 0.001***",
    paste0("p = ", signif(p, 3), stars)
  )
}

## a
ed1a <- ggplot(per_genome, aes(x = asm_total_len_bp / 1e6, y = asm_scaffold_N50_Mb, color = kingdom)) +
  geom_point(alpha = 0.6, size = 1) +
  scale_color_manual(values = kingdom_cols) +
  labs(x = "Assembly size (Mb)", y = "Scaffold N50 (Mb)") +
  theme_nature() +
  theme(
    legend.position = c(0.98, 0.80),
    legend.justification = c("right", "top"),
    legend.background = element_rect(fill = scales::alpha("white", 0.7), colour = "grey"),
    legend.key.size = grid::unit(3, "mm"),
    legend.text = element_text(size = 6),
    legend.spacing.y = grid::unit(1, "mm"),
    legend.title = element_blank()
  )

## b
stat_busco <- per_genome %>%
  rstatix::wilcox_test(BUSCO_completeness ~ kingdom, p.adjust.method = "BH") %>%
  mutate(label = p_label(p.adj))

ymax_busco <- max(per_genome$BUSCO_completeness, na.rm = TRUE)
stat_busco <- stat_busco %>%
  ggpubr::add_y_position(fun = "max", step.increase = 0.05) %>%
  mutate(y.position = ymax_busco * c(1.04, 1.08, 1.12))

ed1b <- ggplot(per_genome, aes(x = kingdom, y = BUSCO_completeness, fill = kingdom, color = kingdom)) +
  geom_jitter(width = 0.15, alpha = 0.6, size = 1.0) +
  geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.75) +
  scale_fill_manual(values = kingdom_cols) +
  scale_color_manual(values = kingdom_cols) +
  labs(x = NULL, y = "BUSCO completeness (%)") +
  stat_pvalue_manual(
    stat_busco,
    label = "label",
    y.position = "y.position",
    tip.length = 0.01,
    bracket.size = 0.3,
    size = 2.0,
    inherit.aes = FALSE
  ) +
  theme_nature()

## c
stat_te <- per_genome %>%
  rstatix::wilcox_test(TE_total_pct ~ kingdom, p.adjust.method = "BH") %>%
  mutate(label = p_label(p.adj))

ymax_te <- max(per_genome$TE_total_pct, na.rm = TRUE)
stat_te <- stat_te %>%
  ggpubr::add_y_position(fun = "max", step.increase = 0.05) %>%
  mutate(y.position = ymax_te * c(1.10, 1.18, 1.26))

ed1c <- ggplot(per_genome, aes(x = kingdom, y = TE_total_pct, fill = kingdom, color = kingdom)) +
  geom_jitter(width = 0.15, alpha = 0.6, size = 1.0) +
  geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.75) +
  scale_fill_manual(values = kingdom_cols) +
  scale_color_manual(values = kingdom_cols) +
  labs(x = NULL, y = "TE proportion (%)") +
  stat_pvalue_manual(
    stat_te,
    label = "label",
    y.position = "y.position",
    tip.length = 0.01,
    bracket.size = 0.3,
    size = 2.0,
    inherit.aes = FALSE
  ) +
  theme_nature()

## d
ed1d <- ggplot(per_genome, aes(x = gene_density, y = kingdom, fill = kingdom)) +
  ggridges::geom_density_ridges(alpha = 0.8) +
  scale_fill_manual(values = kingdom_cols) +
  labs(x = "Gene density (genome)", y = NULL) +
  theme_nature() +
  theme(axis.text.y = element_text(angle = 90, hjust = 0.5, vjust = 0.5))

## e
ed1e <- ggplot(per_genome, aes(x = pi_genome_wide, color = kingdom)) +
  geom_density(linewidth = 0.8) +
  scale_color_manual(values = kingdom_cols) +
  labs(x = "Genome-wide nucleotide diversity (π)", y = "Density") +
  theme_nature()

## f
ed1f <- ggplot(per_genome, aes(x = indel_rate, color = kingdom)) +
  geom_density(linewidth = 0.8) +
  scale_color_manual(values = kingdom_cols) +
  labs(x = "Genome-wide indel rate", y = "Density") +
  guides(color = guide_legend(title = NULL)) +
  theme_nature() +
  theme(
    legend.position = c(0.55, 0.90),
    legend.justification = c(0, 1),
    legend.direction = "vertical",
    legend.text = element_text(size = 7),
    legend.key = element_rect(fill = "transparent", colour = NA),
    legend.key.height = grid::unit(3, "mm"),
    legend.key.width = grid::unit(3, "mm")
  )

Extended_Data_Fig1 <- (ed1a | ed1b | ed1c) /
  (ed1d | ed1e | ed1f) +
  plot_annotation(tag_levels = "a")

save_pdf(Extended_Data_Fig1, "Extended_Data_Fig1.pdf", 180, 120)
save_png(Extended_Data_Fig1, "Extended_Data_Fig1.png", 180, 120, dpi = 300)

message("Done: Extended_Data_Fig1.pdf / Extended_Data_Fig1.png")
