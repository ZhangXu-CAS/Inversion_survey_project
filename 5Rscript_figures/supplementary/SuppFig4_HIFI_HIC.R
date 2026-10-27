#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
})
OUT_DIR <- file.path(PROJ_ROOT, 'outputs', 'supplementary')
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

set.seed(20260925)
NBOOT <- 5000L

pg <- read_csv("source_audit/per_genome.csv", show_col_types = FALSE)
pv <- read_tsv("source_audit/sample_provenance_resolved.tsv", show_col_types = FALSE)
jr <- read_tsv("source_audit/junction_observed_expected.tsv", show_col_types = FALSE)
vm <- read_tsv("source_audit/validation_model_results.tsv", show_col_types = FALSE)
ce <- read_tsv("source_audit/core_effect_sizes.tsv", show_col_types = FALSE)

kingdom_cols <- c(Fungi = "#3274A1", Metazoa = "#E1812C", Viridiplantae = "#00A087")
prov_cols <- c(Same = "#595959", Different = "#E1812C", `Multi-source` = "#3274A1", Unresolved = "#BDBDBD")

theme_nature <- theme_classic(base_size = 8) +
  theme(
    axis.line = element_line(linewidth = 0.45, colour = "black"),
    axis.ticks = element_line(linewidth = 0.40, colour = "black"),
    axis.title = element_text(size = 8.4, margin = margin(t = 3, r = 3, b = 3, l = 3)),
    axis.text = element_text(size = 7.4, colour = "black"),
    legend.title = element_blank(),
    legend.text = element_text(size = 7),
    legend.key.height = unit(3.2, "mm"),
    legend.key.width = unit(4.5, "mm"),
    strip.background = element_blank(),
    strip.text = element_text(size = 8.1, colour = "black"),
    plot.margin = margin(5, 8, 5, 7),
    plot.tag = element_text(size = 10.5, face = "bold"),
    plot.tag.position = c(0.01, 0.99)
  )

prov <- pg %>%
  select(sample_id, kingdom, major_clade, pi_genome_wide, inv_count_per_Gb, inv_len_median_bp) %>%
  inner_join(pv %>% select(sample_id, final_provenance), by = "sample_id") %>%
  mutate(
    source = case_when(
      final_provenance == "SAME_SOURCE_INDIVIDUAL" ~ "Same",
      final_provenance %in% c("DIFFERENT_SOURCE_INDIVIDUAL", "SAME_GENOTYPE_DIFFERENT_INDIVIDUAL") ~ "Different",
      final_provenance %in% c("MULTI_HIFI_SOURCE", "MULTI_HIC_SOURCE", "MULTI_BOTH_SOURCE") ~ "Multi-source",
      TRUE ~ "Unresolved"
    ),
    source = factor(source, levels = c("Same", "Different", "Multi-source", "Unresolved"))
  )

a.dat <- bind_rows(
  prov %>% count(group = kingdom, source),
  prov %>% filter(kingdom == "Metazoa", major_clade == "Insects") %>% count(group = "Insects", source)
) %>%
  filter(group %in% c("Fungi", "Metazoa", "Insects", "Viridiplantae")) %>%
  group_by(group) %>%
  mutate(percent = 100 * n / sum(n)) %>%
  ungroup() %>%
  mutate(group = factor(group, levels = rev(c("Fungi", "Metazoa", "Insects", "Viridiplantae"))))
a.dat$source <- factor(a.dat$source, levels = c("Unresolved", "Same", "Multi-source", "Different"))
a <- ggplot(a.dat, aes(percent, group, fill = source)) +
  geom_col(width = 0.68) +
  scale_fill_manual(values = prov_cols) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.01))) +
  labs(x = "Genomes (%)", y = NULL) +
  theme_nature +
  theme(legend.position = c(0.3,1.1), 
        legend.text = element_text(size = 5),
        legend.key.size = unit(2.5, "mm"),
        legend.direction = "horizontal",
        legend.background = element_blank(),   
        legend.key = element_blank()
      )

boot_diff <- function(dat, value) {
  d <- dat %>% filter(source %in% c("Same", "Different"), is.finite(.data[[value]]))
  s <- d %>% filter(source == "Same") %>% pull(all_of(value))
  z <- d %>% filter(source == "Different") %>% pull(all_of(value))
  if (length(s) < 2 || length(z) < 2) return(tibble(effect = NA, lo = NA, hi = NA, n_diff = length(z), n_same = length(s)))
  b <- replicate(NBOOT, median(sample(z, length(z), TRUE)) - median(sample(s, length(s), TRUE)))
  tibble(effect = median(z) - median(s), lo = quantile(b, 0.05), hi = quantile(b, 0.95), n_diff = length(z), n_same = length(s))
}

effect_input <- prov %>%
  mutate(
 #   log10_pi = if_else(pi_genome_wide > 0, log10(pi_genome_wide), NA_real_),
    log10_count = if_else(inv_count_per_Gb > 0, log10(inv_count_per_Gb), NA_real_),
    log10_length = if_else(inv_len_median_bp > 0, log10(inv_len_median_bp), NA_real_)
  )

groups <- list(
  Metazoa = effect_input %>% filter(kingdom == "Metazoa"),
  Insects = effect_input %>% filter(kingdom == "Metazoa", major_clade == "Insects"),
  Viridiplantae = effect_input %>% filter(kingdom == "Viridiplantae")
)

metrics <- c(
 # log10_pi = expression(pi[genome]),
  log10_count = "Inversions/Gb",
  log10_length = "Median inversion\nlength"
)

b.dat <- imap_dfr(groups, function(dat, grp) {
  map_dfr(names(metrics), function(m) {
    boot_diff(dat, m) %>% mutate(group = grp, metric = m)
  })
}) %>%
  mutate(
    group = factor(group, levels = rev(c("Metazoa", "Insects", "Viridiplantae"))),
    colour_group = if_else(group == "Viridiplantae", "Viridiplantae", "Metazoa"),
    metric = factor(metric, levels = names(metrics), labels = c( "Inversions/Gb", "Median inversion length")),
    n = paste0(n_diff, "/", n_same)
  )

write_tsv(b.dat, "source_effects_three_metrics.tsv")

b <- ggplot(b.dat, aes(effect, group, colour = colour_group)) +
  geom_vline(xintercept = 0, linetype = 2, linewidth = 0.40) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0, linewidth = 0.70) +
  geom_point(size = 2.15) +
  facet_wrap(~metric, nrow = 1, scales = "free_x") +
  geom_text(aes(x = hi, label = n), hjust = -0.12, colour = "black", size = 2.15) +
  scale_colour_manual(values = kingdom_cols) +
  labs(x = "Different − same", y = NULL) +
  theme_nature +
  theme(legend.position = "none",panel.spacing.x = unit(1.5, "lines"), 
        plot.margin = margin(t = 5.5, r = 20, b = 5.5, l = 5.5, unit = "pt") ) +
  coord_cartesian(clip = "off")

jg <- jr %>%
  left_join(pg %>% select(sample_id, kingdom), by = "sample_id") %>%
  group_by(sample_id, kingdom, threshold_bp) %>%
  summarise(
    O = mean(O_d, na.rm = TRUE),
    E = mean(E_d, na.rm = TRUE),
    null = mean(null_mean, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(across(c(O, E, null), ~ifelse(is.nan(.x), NA_real_, .x))) %>%
  filter(is.finite(O), is.finite(E), E > 0, is.finite(null))

boot_ratio <- function(dat) {
  n <- nrow(dat)
  obs <- mean(dat$O) / mean(dat$E)
  b <- replicate(NBOOT, {
    i <- sample.int(n, n, TRUE)
    mean(dat$O[i]) / mean(dat$E[i])
  })
  tibble(ratio = obs, lo = quantile(b, 0.025), hi = quantile(b, 0.975), N = n)
}

c.dat <- jg %>%
  filter(kingdom %in% c("Metazoa", "Viridiplantae"), threshold_bp %in% c(10000, 50000)) %>%
  group_by(kingdom, threshold_bp) %>%
  group_modify(~boot_ratio(.x)) %>%
  ungroup() %>%
  mutate(
    threshold = factor(threshold_bp, levels = c(10000, 50000), labels = c("10 kb", "50 kb")),
    y = log2(ratio), yl = log2(lo), yh = log2(hi)
  )

write_tsv(c.dat, "junction_kingdom_summary_for_figure.tsv")

c <- ggplot(c.dat, aes(threshold, y, colour = kingdom, group = kingdom)) +
  geom_hline(yintercept = 0, linetype = 2, linewidth = 0.40) +
  geom_line(linewidth = 0.65, position = position_dodge(0.14)) +
  geom_errorbar(aes(ymin = yl, ymax = yh), width = 0.07, linewidth = 0.60, position = position_dodge(0.14)) +
  geom_point(size = 2.10, position = position_dodge(0.14)) +
  scale_colour_manual(values = kingdom_cols) +
  labs(x = "Breakpoint distance", y = expression(log[2]*"(observed / expected)")) +
  theme_nature +
  theme(legend.position = "top", 
        axis.title.y = element_text( margin = margin(t = 0, r = -1, b = 0, l = 0)), 
        legend.box.margin = margin(0, 0, -2, 0))

d.dat <- vm %>%
  filter(threshold_kb %in% c(10, 50), term %in% c("different", "near", "different:near")) %>%
  mutate(
    term = recode(
      term,
      different = "Different\nsource",
      near = "Near scaffold\njoin",
      `different:near` = "Different source\n× near join"
    ),
    term = factor(term, levels = rev(c("Different\nsource", "Near scaffold\njoin", "Different source\n× near join"))),
    threshold = factor(threshold_kb, levels = c(10, 50), labels = c("10 kb", "50 kb"))
  )

d <- ggplot(d.dat, aes(odds_ratio, term, shape = threshold, group = threshold)) +
  geom_vline(xintercept = 1, linetype = 2, linewidth = 0.40) +
  geom_errorbarh(aes(xmin = OR_CI_low, xmax = OR_CI_high), height = 0, linewidth = 0.62, position = position_dodge(0.30)) +
  geom_point(size = 2.05, position = position_dodge(0.30)) +
  scale_x_log10(breaks = c(0.2, 0.5, 1, 2), labels = c("0.2", "0.5", "1", "2")) +
  labs(x = "Odds ratio", y = NULL) +
  theme_nature +
  theme(legend.position = "top", legend.box.margin = margin(0, 0, -2, 0))

set_order <- c("S0_FULL", "S1_SAME_SOURCE", "S7_SAME_STRICT_INTERNAL_50KB", "S8_NO_MULTI_HIFI")
set_labels <- c(
  S0_FULL = "Full",
  S1_SAME_SOURCE = "Same source",
  S7_SAME_STRICT_INTERNAL_50KB = "Same\n(internal >50 kb)",
  S8_NO_MULTI_HIFI = "No multi-HiFi"
)

e.dat <- ce %>%
  filter(dataset %in% set_order, metric %in% c("log10_inversion_per_Gb", "log10_median_inversion_length")) %>%
  mutate(
    filter_set = factor(unname(set_labels[dataset]), levels = rev(unname(set_labels[set_order]))),
    metric = recode(
      metric,
      log10_inversion_per_Gb = "Inversion number",
      log10_median_inversion_length = "Median inversion length"
    ),
    metric = factor(metric, levels = c("Inversion density", "Median inversion length"))
  )

e <- ggplot(e.dat, aes(effect, filter_set, colour = metric)) +
  geom_vline(xintercept = 0, linetype = 2, linewidth = 0.40) +
  geom_errorbarh(aes(xmin = CI_low, xmax = CI_high), height = 0, linewidth = 0.72) +
  geom_point(size = 2.20) +
  facet_wrap(~metric, nrow = 1) +
  scale_colour_manual(values = c(
    "Inversion number" = kingdom_cols["Metazoa"],
    "Median inversion length" = kingdom_cols["Viridiplantae"]
  )) +
  scale_x_continuous(limits = c(0, 0.65), breaks = c(0, 0.2, 0.4, 0.6), expand = expansion(mult = c(0.02, 0.04))) +
  labs(x = "Cross-kingdom effect", y = NULL) +
  theme_nature +
  theme(legend.position = "none")

design <- "
AAAABBBBBB
CCCCCDDDDD
EEEEEEEEEE
"


fig <- a + b + c + d + e +
  plot_layout(design = design) +
  plot_annotation(tag_levels = "a")

ggsave("outputs/supplementary/SuppFig4_HiFi_HiC.pdf", fig, width = 180, height = 170, units = "mm", device = cairo_pdf)
ggsave("outputs/supplementary/SuppFig4_HiFi_HiC.png", fig, width = 180, height = 170, units = "mm", dpi = 600)
