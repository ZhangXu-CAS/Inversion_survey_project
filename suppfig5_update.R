############################################################
## Supplementary Figure Robustness to genome-wide heterozygosity filter
##
## Rationale:
##   Very low genome-wide nucleotide diversity (pi_genome_wide) can indicate
##   near-isogenic/selfing lineages (especially in plants), which may affect
##   inversion detection and/or the interpretation of inversion "age" proxies.
############################################################

  library(tidyverse)
  library(scales)
  library(patchwork)
## =========================
## 0) User-editable settings
## =========================
IN_PER_GENOME    <- "per_genome.csv"
IN_PER_INVERSION <- "per_inversion.csv"
OUT_DIR          <- "SuppFig5-6_pi_filter_outputs"


## Age threshold for "old" inversions (π_inversion / π_genome)
AGE_THRESHOLD <- 1.0

## Downsampling for scatter (per kingdom)
MAX_POINTS_PER_KINGDOM <- 10000L
SEED <- 123

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

save_pdf <- function(plot, filename, width_mm, height_mm) {
  ggsave(
    filename = file.path(OUT_DIR, filename),
    plot     = plot,
    width    = width_mm,
    height   = height_mm,
    units    = "mm",
    useDingbats = FALSE
  )
}

save_png <- function(plot, filename, width_mm, height_mm, dpi = 300) {
  ggsave(
    filename = file.path(OUT_DIR, filename),
    plot     = plot,
    width    = width_mm,
    height   = height_mm,
    units    = "mm",
    dpi      = dpi
  )
}
## =========================
##  Style
## =========================
kingdom_cols <- c(
  "Fungi"         = "#3274A1",   # kept for global consistency (unused here)
  "Metazoa"       = "#E1812C",
  "Viridiplantae" = "#00A087"
)

theme_nature <- function(base_size = 8) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid      = element_blank(),
      axis.text       = element_text(color = "black"),
      axis.title      = element_text(color = "black"),
      legend.title    = element_blank(),
      panel.border    = element_rect(color = "black", linewidth = 0.5),
      axis.ticks      = element_line(color = "black"),
      plot.title      = element_blank(),
      plot.subtitle   = element_blank(),
      plot.margin     = margin(1.5, 1.5, 1.5, 1.5, unit = "mm")
    )
}

## =========================
##  Load and harmonize columns
## =========================
per_inv_raw <- readr::read_csv(IN_PER_INVERSION, show_col_types = FALSE)
## ---- Load per-genome table and apply pi_genome_wide filter ----
per_genome_raw <- readr::read_csv(IN_PER_GENOME, show_col_types = FALSE)

## Identify a shared sample/genome identifier column between the two tables
id_candidates <- c("sample_id", "genome_id", "genome", "assembly_id", "assembly", "accession", "acc", "id")
id_col <- intersect(id_candidates, intersect(names(per_genome_raw), names(per_inv_raw)))[1]
if (is.na(id_col)) {
  stop(
    "Cannot find a shared ID column between per_genome.csv and per_inversion.csv.\n",
    "Tried: ", paste(id_candidates, collapse = ", "), "\n",
    "per_genome columns: ", paste(names(per_genome_raw), collapse = ", "), "\n",
    "per_inversion columns: ", paste(names(per_inv_raw), collapse = ", ")
  )
}

## Identify pi_genome_wide column (allow a few common variants)
pi_candidates <- c("pi_genome_wide", "pi_genomewide", "pi_genome", "pi_genome_all", "pi_genomewide_mean")
pi_col <- intersect(pi_candidates, names(per_genome_raw))[1]
if (is.na(pi_col)) {
  stop(
    "Cannot find genome-wide pi column in per_genome.csv.\n",
    "Expected one of: ", paste(pi_candidates, collapse = ", "), "\n",
    "per_genome columns: ", paste(names(per_genome_raw), collapse = ", ")
  )
}

## Filter threshold
PI_GENOME_MIN4 <- 1e-4
PI_GENOME_MIN3 <- 1e-3

###PI_GENOME_MIN4 <- 1e-4
per_genome_filt4 <- per_genome_raw %>%
  mutate(
    pi_genome_wide = as.numeric(.data[[pi_col]])
  ) %>%
  filter(!is.na(.data[[id_col]])) %>%
  filter(!is.na(pi_genome_wide))%>%
  filter(pi_genome_wide >= PI_GENOME_MIN4)

keep_ids4 <- per_genome_filt4 %>%
  filter(pi_genome_wide >= PI_GENOME_MIN4) %>%
  pull(.data[[id_col]]) %>%
  unique()
message("  per_genome: ", nrow(per_genome_raw), " rows -> ", length(keep_ids4), " genomes kept after pi> 1e-4")

## Save filtered tables for transparency
readr::write_csv(
  per_genome_filt4 %>% filter(.data[[id_col]] %in% keep_ids4),
  file.path(OUT_DIR, "per_genome_pi_ge_1e4.csv")
)
readr::write_csv(per_inv_raw, file.path(OUT_DIR, "per_inversion_pi_ge_1e4.csv"))
#######
###PI_GENOME_MIN4 <- 1e-3
per_genome_filt3 <- per_genome_raw %>%
  mutate(
    pi_genome_wide = as.numeric(.data[[pi_col]])
  ) %>%
  filter(!is.na(.data[[id_col]])) %>%
  filter(!is.na(pi_genome_wide))%>%
  filter(pi_genome_wide >= PI_GENOME_MIN3)

keep_ids3 <- per_genome_filt3 %>%
  filter(pi_genome_wide >= PI_GENOME_MIN3) %>%
  pull(.data[[id_col]]) %>%
  unique()
message("  per_genome: ", nrow(per_genome_raw), " rows -> ", length(keep_ids3), " genomes kept after pi> 1e-3")

## Save filtered tables for transparency
readr::write_csv(
  per_genome_filt3 %>% filter(.data[[id_col]] %in% keep_ids3),
  file.path(OUT_DIR, "per_genome_pi_ge_1e3.csv")
)

## Helper: p-value to stars
p_to_stars <- function(p) {
  dplyr::case_when(
    is.na(p)      ~ "",
    p < 0.0001    ~ "****",
    p < 0.001     ~ "***",
    p < 0.01      ~ "**",
    p < 0.05      ~ "*",
    TRUE          ~ "ns"
  )
}
####update liner regrssion with pi (filtered pi >1e-4 ) --03-03-2026
#-------------------------
# Robust numeric coercion 
#-------------------------
to_num <- function(x){
  suppressWarnings(readr::parse_number(as.character(x)))
}

#-------------------------
# BH + p label:  p<0.001 -> "p < 0.001***"
# else -> "p = X" + stars (**/*; >=0.05 no star)
#-------------------------
p_label_vec <- function(p){
  stars <- dplyr::case_when(
    is.na(p) ~ "",
    p < 1e-2 ~ "**",
    p < 5e-2 ~ "*",
    TRUE     ~ ""
  )
  dplyr::case_when(
    is.na(p) ~ "p = NA",
    p < 1e-3 ~ "p < 0.001***",
    TRUE     ~ paste0("p = ", signif(p, 3), stars)
  )
}

#-------------------------
# Build per-kingdom legend labels: "Kingdom \n R²=..., p=..."
# stats: Pearson p (or lm slope p if you switch)
#-------------------------
make_legend_labels_r2p <- function(df, x, y, method = c("pearson","lm_slope"), p_adjust = TRUE){
  method <- match.arg(method)
  
  stat <- df %>%
    group_by(kingdom) %>%
    summarise(
      n = sum(is.finite(.data[[x]]) & is.finite(.data[[y]])),
      sd_x = sd(.data[[x]], na.rm = TRUE),
      sd_y = sd(.data[[y]], na.rm = TRUE),
      r = ifelse(n >= 3 && sd_x > 0 && sd_y > 0,
                 cor(.data[[x]], .data[[y]], use = "complete.obs", method = "pearson"),
                 NA_real_),
      p = ifelse(n >= 3 && sd_x > 0 && sd_y > 0,
                 {
                   if (method == "pearson") {
                     cor.test(.data[[x]], .data[[y]], method = "pearson")$p.value
                   } else {
                     # lm slope p-value
                     summary(lm(.data[[y]] ~ .data[[x]]))$coefficients[2, 4]
                   }
                 },
                 NA_real_),
      .groups = "drop"
    )
  
  # BH adjust across kingdoms (non-NA only)
  stat$p_adj <- NA_real_
  idx <- which(!is.na(stat$p))
  if (p_adjust && length(idx) > 0) stat$p_adj[idx] <- p.adjust(stat$p[idx], method = "BH")
  if (!p_adjust) stat$p_adj <- stat$p
  
  stat <- stat %>%
    mutate(
      r2 = r^2,
      leg = ifelse(
        is.na(p_adj) | is.na(r2),
        paste0(as.character(kingdom), "\n", "n=", n, " (insufficient/constant)"),
        paste0(as.character(kingdom), "\nR²=", signif(r2, 3), ", ", gsub("^p = ", "p=", p_label_vec(p_adj)))
      )
    )
  
  setNames(stat$leg, as.character(stat$kingdom))
}

#-------------------------
corr_panel_kingdom <- function(df,
                               xvar, yvar,
                               xlab, ylab,
                               logx = FALSE, logy = FALSE,
                               method_stats = "pearson",    
                               legend_pos = c(0.70, 0.20),
                               legend_text_size = 5,
                               point_alpha = 0.5,
                               point_size = 1.0,
                               line_width = 0.8){
  
  
  d <- df %>%
    mutate(
      kingdom = factor(kingdom, levels = c("Fungi","Metazoa","Viridiplantae")),
      .x = to_num(.data[[xvar]]),
      .y = to_num(.data[[yvar]])
    )
  
  # safe log
  eps <- 1e-6
  if (logx) d <- d %>% mutate(.x = log10(pmax(.x, eps)))
  if (logy) d <- d %>% mutate(.y = log10(pmax(.y, eps)))
  
  d <- d %>% filter(!is.na(kingdom), is.finite(.x), is.finite(.y))
  
  legend_labels <- make_legend_labels_r2p(d, ".x", ".y", method = method_stats, p_adjust = TRUE)
  
  ggplot(d, aes(x = .x, y = .y, color = kingdom)) +
    geom_point(alpha = point_alpha, size = point_size) +
    geom_smooth(method = "lm", se = FALSE, linewidth = line_width) +
    scale_color_manual(
      values = kingdom_cols,
      labels = legend_labels,
      guide  = guide_legend(title = NULL)
    ) +
    labs(x = xlab, y = ylab) +
    theme_nature() +
    theme(
      legend.position   = legend_pos,
      legend.text       = element_text(size = legend_text_size),
      legend.key.height = unit(12, "pt"),
      legend.background = element_rect(fill = "transparent", color = NA)
    )
}

# 1 rows pi > 1e-4
S0 <- corr_panel_kingdom(
  per_genome_filt4,
  xvar = "pi_genome_wide",
  yvar = "inv_count_per_Gb",
  xlab = "Genome-wide nucleotide diversity (π)",
  ylab = expression(paste("Inversion count per Gb (log"[10],")")),
  logy = TRUE,              
  method_stats = "pearson",
  legend_pos = c(0.70, 0.20),
  legend_text_size = 5
)

#  inversion counts vs π
S1 <- corr_panel_kingdom(
  per_genome_filt4,
  xvar = "pi_genome_wide",
  yvar = "inv_count_total",
  xlab = "Genome-wide nucleotide diversity (π)",
  ylab = expression(paste("Total inversion count (log"[10],")")),
  logy = TRUE,
  legend_pos = c(0.70, 0.20),
  legend_text_size = 5
)

# median inversion length vs π（median_log10_inv_length）
S2 <- corr_panel_kingdom(
  per_genome_filt4,
  xvar = "pi_genome_wide",
  yvar = "median_log10_inv_length",
  xlab = "Genome-wide nucleotide diversity (π)",
  ylab = expression(paste("Median inversion length (log"[10], " bp)")),
  logy = !("median_log10_inv_length" %in% names(per_genome_filt4)),
  legend_pos = c(0.70, 0.85),
  legend_text_size = 5
)

# 2 rows pi > 1e-3
S3 <- corr_panel_kingdom(
  per_genome_filt3,
  xvar = "pi_genome_wide",
  yvar = "inv_count_per_Gb",
  xlab = "Genome-wide nucleotide diversity (π)",
  ylab = expression(paste("Inversion count per Gb (log"[10],")")),
  logy = TRUE,              
  method_stats = "pearson",
  legend_pos = c(0.70, 0.20),
  legend_text_size = 5
)

#  inversion counts vs π
S4 <- corr_panel_kingdom(
  per_genome_filt3,
  xvar = "pi_genome_wide",
  yvar = "inv_count_total",
  xlab = "Genome-wide nucleotide diversity (π)",
  ylab = expression(paste("Total inversion count (log"[10],")")),
  logy = TRUE,
  legend_pos = c(0.70, 0.20),
  legend_text_size = 5
)

#  median inversion length vs π（median_log10_inv_length）
S5 <- corr_panel_kingdom(
  per_genome_filt3,
  xvar = "pi_genome_wide",
  yvar = "median_log10_inv_length",
  xlab = "Genome-wide nucleotide diversity (π)",
  ylab = expression(paste("Median inversion length (log"[10], " bp)")),
  logy = !("median_log10_inv_length" %in% names(per_genome_filt3)),
  legend_pos = c(0.70, 0.85),
  legend_text_size = 5
)

SuppFig5 <- ( S0 | S1 | S2 ) / (S3 | S4 | S5) +
  plot_annotation(tag_levels = "a",
                  theme = theme(plot.tag = element_text(face = "bold", size = 8)))


## Main exports (Nature-style square-ish multi-panel)
save_pdf(SuppFig5, "SuppFig5.pdf", 180, 120)
save_png(SuppFig5, "SuppFig5.png", 180, 120, dpi = 300)


####SuppFig6 rpeat fig 3
per_inv_raw <- per_inv_raw %>%
  filter(.data[[id_col]] %in% keep_ids4)
message("  per_inversion (after):  ", nrow(per_inv_raw), " rows")

per_inv <- per_inv_raw %>%
  filter(
    kingdom %in% c("Metazoa", "Viridiplantae"),
    !is.na(inv_len_bp)
  ) %>%
  mutate(
    kingdom    = factor(kingdom, levels = c("Metazoa", "Viridiplantae")),
    inv_len_bp = as.numeric(inv_len_bp),
    log10_len  = log10(inv_len_bp),
    
    ## Age proxy
    pi_ratio = as.numeric(pi_ratio_inv_genome),
    pi_ratio = if_else(pi_ratio <= 0 | !is.finite(pi_ratio), NA_real_, pi_ratio),
    
    ## Age class (relative to threshold)
    age_class = case_when(
      is.na(pi_ratio)           ~ NA_character_,
      pi_ratio <= AGE_THRESHOLD ~ "Young",
      pi_ratio >  AGE_THRESHOLD ~ "Old"
    ),
    age_class = factor(age_class, levels = c("Young", "Old")),
    
    ## Size class (S/M/L)
    size_class = case_when(
      inv_len_bp < 1e4                      ~ "Small",
      inv_len_bp >= 1e4 & inv_len_bp <= 1e6 ~ "Medium",
      inv_len_bp > 1e6                      ~ "Large",
      TRUE                                  ~ NA_character_
    ),
    size_class = factor(size_class, levels = c("Small", "Medium", "Large"))
  )

per_inv_age <- per_inv %>% filter(!is.na(pi_ratio), is.finite(log10_len))

## Reasonable plotting range for length (log10)
log10_min <- max(2, floor(min(per_inv_age$log10_len, na.rm = TRUE)))
log10_max <- min(8, ceiling(max(per_inv_age$log10_len, na.rm = TRUE)))

## ---- A: Length vs age proxy (scatter + linear fit; faceted by size class) ----
## Requirements: size_class already defined in per_inv_age as factor Small/Medium/Large

## ---- A: Length vs log10(age proxy) (scatter + LM; faceted by size class) ----
## Key changes:
##  - y: log10(pi_ratio + eps)
##  - annotate only r; p as stars
##  - x-axis breaks: integers only

set.seed(SEED)

EPS_PI <- 1e-6  # small offset to avoid log10(0); tweak if your pi_ratio can be extremely small

## ---- Length vs age proxy (scatter + LOESS; downsampled; faceted) ----
set.seed(SEED)
dat_A <- per_inv_age 
dat_A_filter <- per_inv_age %>%
  filter(log10_len >= 3, log10_len <= 8) %>%
  group_by(kingdom) %>%
  group_modify(~{
    n_here <- nrow(.x)
    n_take <- min(MAX_POINTS_PER_KINGDOM, n_here)
    if (n_here > n_take) slice_sample(.x, n = n_take) else .x
  }) %>%
  ungroup()

pi99_A <- quantile(dat_A$pi_ratio, 0.99, na.rm = TRUE)

SuppFig6A <- ggplot(dat_A, aes(x = log10_len, y = pi_ratio)) +
  geom_point(aes(color = kingdom), alpha = 0.10, size = 0.3, show.legend = FALSE) +
  geom_smooth(data=dat_A_filter,aes(color = kingdom), method = "loess", se = TRUE, linewidth = 0.8, span = 0.7, show.legend = FALSE) +
  scale_color_manual(values = kingdom_cols) +
  facet_wrap(~ kingdom, nrow = 1) +
  coord_cartesian(xlim = c(2, 8),
                  ylim = c(0, pi99_A)
  ) +
  geom_hline(yintercept = AGE_THRESHOLD, linetype = "dashed", linewidth = 0.35) +
  labs(
    x = expression(paste("Inversion length (log"[10], " bp)")),
    y = expression(pi[inversion] / pi[genome])
  ) +
  theme_nature()

## ---- B: Kingdom-normalized 2D relative density (length × age) ----
dat_len_age <- per_inv_age %>%
  filter(
    is.finite(log10_len),
    is.finite(pi_ratio),
    log10_len >= log10_min,
    log10_len <= log10_max,
    kingdom %in% c("Metazoa", "Viridiplantae")
  ) %>%
  mutate(
    kingdom = factor(kingdom, levels = c("Metazoa", "Viridiplantae"))
  )

## shared y scale across A/B
age_max_plot <- quantile(dat_len_age$pi_ratio, 0.99, na.rm = TRUE)
## 
n_len_bins_B <- 50
n_age_bins_B <- 50

len_breaks_B <- seq(log10_min, log10_max, length.out = n_len_bins_B + 1)
age_breaks_B <- seq(0, age_max_plot, length.out = n_age_bins_B + 1)

dat_B_grid <- dat_len_age %>%
  mutate(
    len_bin = cut(log10_len, breaks = len_breaks_B, include.lowest = TRUE, labels = FALSE),
    age_bin = cut(pi_ratio,  breaks = age_breaks_B, include.lowest = TRUE, labels = FALSE)
  ) %>%
  drop_na(len_bin, age_bin) %>%
  count(kingdom, len_bin, age_bin, name = "n") %>%
  group_by(kingdom) %>%
  mutate(
    prop     = n / sum(n),
    prop_rel = prop / max(prop),
    ## gamma-like rescaling to sharpen hotspots
    prop_plot = sqrt(prop_rel),
    len_mid  = (len_breaks_B[len_bin] + len_breaks_B[len_bin + 1]) / 2,
    age_mid  = (age_breaks_B[age_bin] + age_breaks_B[age_bin + 1]) / 2
  ) %>%
  ungroup()

SuppFig6B <- ggplot(dat_B_grid, aes(x = len_mid, y = age_mid)) +
  geom_tile(aes(fill = prop_plot)) +
  geom_hline(
    yintercept = AGE_THRESHOLD,
    linetype = "dashed",
    linewidth = 0.35,
    colour = "black"
  ) +
  facet_wrap(~ kingdom, nrow = 1) +
  coord_cartesian(
    xlim = c(log10_min, log10_max),
    ylim = c(0, age_max_plot)
  ) +
  scale_fill_gradientn(
    colours = c("#f7fbff", "#deebf7", "#9ecae1", "#3182bd", "#31a354", "#ffff00"),
    values  = rescale(c(0, 0.12, 0.28, 0.5, 0.75, 1)),
    limits  = c(0, 1),
    name    = "Relative density",
    guide   = guide_colorbar(
      barheight = unit(15, "mm"),
      barwidth  = unit(3, "mm"),
      ticks.colour = "black",
      frame.colour = "black"
    )
  ) +
  labs(
    x = expression(paste("Inversion length (log"[10], " bp)")),
    y = expression(pi[inversion] / pi[genome])
  ) +
  theme_nature() +
  theme(
    legend.title     = element_text(size = 5, lineheight = 0.9),
    legend.text      = element_text(size = 5),
    legend.key.size  = unit(3, "mm")
  )

## ---- C: Logistic model P(Old) vs length + 95% CI ----
dat_glm <- per_inv_age %>%
  filter(!is.na(age_class), log10_len >= log10_min, log10_len <= log10_max)

glm_old <- glm(
  I(age_class == "Old") ~ log10_len * kingdom,
  data   = dat_glm,
  family = binomial()
)

newdat <- expand.grid(
  log10_len = seq(log10_min, log10_max, length.out = 200),
  kingdom   = c("Metazoa", "Viridiplantae")
)

pred <- predict(glm_old, newdata = newdat, se.fit = TRUE)

newdat <- newdat %>%
  mutate(
    fit    = pred$fit,
    se     = pred$se.fit,
    p      = plogis(fit),
    p_low  = plogis(fit - 1.96 * se),
    p_high = plogis(fit + 1.96 * se)
  )

SuppFig6C <- ggplot(newdat, aes(x = log10_len, y = p, colour = kingdom, fill = kingdom)) +
  geom_ribbon(aes(ymin = p_low, ymax = p_high), alpha = 0.2, colour = NA, show.legend = TRUE) +
  geom_line(linewidth = 0.9, show.legend = TRUE) +
  scale_colour_manual(values = kingdom_cols) +
  scale_fill_manual(values = kingdom_cols) +
  coord_cartesian(xlim = c(log10_min, log10_max), ylim = c(0, 1)) +
  labs(
    x = expression(paste("Inversion length (log"[10], " bp)")),
    y = "Predicted P(old)"
  ) +
  theme_nature()+
  theme(
    legend.position   =c(0.75, 0.2), legend.background = element_rect(fill = NA),
    legend.key.size = unit(3, "mm"),
    legend.title      = element_text(size = 5, lineheight = 0.9),
    legend.text       = element_text(size = 5),
    legend.margin     = margin(0, 0, 0, 0, unit = "mm"),
    legend.box.margin = margin(0, 0, 0, 0, unit = "mm")
  )

## ---- D: Length-binned fraction of Old inversions ----
bin_width <- 0.5
n_bins <- ceiling((log10_max - log10_min) / bin_width)

dat_D <- per_inv_age %>%
  filter(log10_len >= log10_min, log10_len <= log10_max, !is.na(age_class)) %>%
  mutate(
    bin_index = pmin(pmax(floor((log10_len - log10_min) / bin_width), 0L), n_bins - 1L),
    len_mid   = log10_min + bin_width * (bin_index + 0.5)
  ) %>%
  group_by(kingdom, len_mid) %>%
  summarise(
    n_total  = n(),
    n_old    = sum(age_class == "Old"),
    frac_old = n_old / n_total,
    .groups  = "drop"
  )

SuppFig6D <- ggplot(dat_D, aes(x = len_mid, y = frac_old, color = kingdom)) +
  geom_line(linewidth = 0.9, show.legend = TRUE) +
  geom_point(size = 1.4, show.legend = TRUE) +
  scale_color_manual(values = kingdom_cols) +
  coord_cartesian(xlim = c(log10_min, log10_max), ylim = c(0, 1)) +
  labs(
    x = expression(paste("Inversion length (log"[10], " bp, binned)")),
    y = "Fraction of old inversions"
  ) +
  theme_nature()+
  theme(
    legend.position   =c(0.75, 0.2), legend.background = element_rect(fill = NA),
    legend.key.size = unit(3, "mm"),
    legend.title      = element_text(size = 5, lineheight = 0.9),
    legend.text       = element_text(size = 5),
    legend.margin     = margin(0, 0, 0, 0, unit = "mm"),
    legend.box.margin = margin(0, 0, 0, 0, unit = "mm")
  )

## ---- E: Age proxy by size class (scatter + box) + within-class significance ----
dat_E <- per_inv_age %>%
  filter(!is.na(size_class), log10_len >= log10_min, log10_len <= log10_max)

pi99_E <- quantile(dat_E$pi_ratio, 0.99, na.rm = TRUE)

## within each size class: plants vs animals (Wilcoxon rank-sum)
stat_E <- dat_E %>%
  group_by(size_class) %>%
  summarise(
    p_value = tryCatch(
      wilcox.test(pi_ratio ~ kingdom)$p.value,
      error = function(e) NA_real_
    ),
    .groups = "drop"
  ) %>%
  mutate(
    sig = p_to_stars(p_value)
  )

## y-position for significance labels (kept inside axis limits)
## compute per size class to avoid overlap
pos_E <- dat_E %>%
  group_by(size_class) %>%
  summarise(y = quantile(pi_ratio, 0.98, na.rm = TRUE), .groups = "drop") %>%
  mutate(y = pmin(y + 0.06 * pi99_E, pi99_E * 0.98))

stat_E <- stat_E %>% left_join(pos_E, by = "size_class")

SuppFig6E <- ggplot(dat_E, aes(x = size_class, y = pi_ratio,  fill = kingdom)) +
  geom_point(aes(color = kingdom),
             position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.7),
             alpha = 0.25,
             size = 0.45,
             show.legend = TRUE
  ) +
  geom_boxplot(
    position = position_dodge(width = 0.7),
    width = 0.45,
    outlier.shape = NA,
    linewidth = 0.35,
    show.legend = TRUE
  ) +
  geom_hline(yintercept = AGE_THRESHOLD, linetype = "dashed", linewidth = 0.35) +
  geom_text(
    data = stat_E,
    aes(x = size_class, y = y, label = sig),
    inherit.aes = FALSE,
    size = 3.2,
    vjust = 0
  ) +
  scale_color_manual(values = kingdom_cols) +
  scale_fill_manual(values = kingdom_cols) +
  coord_cartesian(ylim = c(0, pi99_E)) +
  labs(
    x = "Inversion length class",
    y = expression(pi[inversion] / pi[genome])
  ) +
  theme_nature()+
  theme( plot.margin     = margin(1.5, 0, 1.5, 1.5, unit = "mm"))+
  theme(
    legend.position   ="right", legend.key.size = unit(3, "mm"),
    legend.title      = element_text(size = 5, lineheight = 0.9),
    legend.text       = element_text(size = 5),
    legend.margin     = margin(0, 0, 0, -0.75, unit = "mm"),
    legend.box.margin = margin(0, 0, 0, -0.75, unit = "mm")
  )
## ---- F: ECDF of age proxy (inversions ≥10 kb) ----
dat_F <- per_inv_age %>%
  filter(
    inv_len_bp >= 1e4, 
    log10_len >= log10_min, log10_len <= log10_max)

pi95_F <- quantile(dat_F$pi_ratio, 0.95, na.rm = TRUE)

SuppFig6F <- ggplot(dat_F, aes(x = pi_ratio, colour = kingdom)) +
  stat_ecdf(aes(y = after_stat(1 - y)), linewidth = 0.7, show.legend = TRUE) +
  scale_colour_manual(values = kingdom_cols) +
  coord_cartesian(xlim = c(0, pi95_F)) +
  labs(
    x = expression(pi[inversion] / pi[genome]),
    y = expression(paste("Fraction of inversions ≥ ", pi[inversion]/pi[genome]))
  ) +
  theme_nature()+
  theme(
    legend.position   =c(0.75,0.8),axis.title.y = element_text(size = 7),
    legend.background = element_rect(fill = NA),
    legend.key.size = unit(3, "mm"),
    legend.title      = element_text(size = 5, lineheight = 0.9),
    legend.text       = element_text(size = 5),
    legend.margin     = margin(0, 0, 0, 0, unit = "mm"),
    legend.box.margin = margin(0, 0, 0, 0, unit = "mm")
  )

## ---- G: Sensitivity of Old fraction to threshold (Large only; ≥1 Mb) ----
threshold_grid <- seq(0.3, 1.5, by = 0.1)

dat_large <- per_inv_age %>%
  filter(
    size_class == "Large", 
    log10_len >= log10_min, log10_len <= log10_max)

sens_res <- purrr::map_dfr(threshold_grid, function(thr) {
  dat_large %>%
    mutate(is_old_thr = pi_ratio > thr) %>%
    group_by(kingdom) %>%
    summarise(
      n_total  = n(),
      frac_old = if_else(n_total > 0, mean(is_old_thr, na.rm = TRUE), NA_real_),
      .groups  = "drop"
    ) %>%
    mutate(age_threshold = thr)
})

SuppFig6G <- ggplot(sens_res, aes(x = age_threshold, y = frac_old, colour = kingdom)) +
  geom_line(linewidth = 0.9, show.legend = TRUE) +
  geom_point(size = 1.2, show.legend = TRUE) +
  scale_colour_manual(values = kingdom_cols) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    x = "Different thresholds for old",
    y = "Fraction of old inversions (≥1 Mb)"
  ) +
  theme_nature()+
  theme(
    legend.position   =c(0.3,0.2), axis.title.y = element_text(size = 7),
    legend.key.size = unit(3, "mm"),
    legend.background = element_rect(fill = NA),
    legend.title      = element_text(size = 5, lineheight = 0.9),
    legend.text       = element_text(size = 5),
    legend.margin     = margin(0, 0, 0, 0, unit = "mm"),
    legend.box.margin = margin(0, 0, 0, 0, unit = "mm")
  )

## ---- H: Bootstrap-resampled fraction of Old inversions by size class ----
dat_H <- per_inv_age %>%
  filter(!is.na(size_class), !is.na(age_class))

n_boot <- 1000

resample_one_size <- function(size_label, n_boot = 1000) {
  dat_sc <- dat_H %>% filter(size_class == size_label)
  tab_sc <- table(dat_sc$kingdom)
  if (length(tab_sc) < 2) return(NULL)
  
  n_min <- min(tab_sc)
  
  purrr::map_dfr(seq_len(n_boot), function(b) {
    dat_b <- dat_sc %>%
      group_by(kingdom) %>%
      slice_sample(n = n_min, replace = FALSE) %>%
      ungroup()
    
    dat_b %>%
      group_by(kingdom) %>%
      summarise(frac_old = mean(age_class == "Old"), .groups = "drop") %>%
      mutate(size_class = size_label, boot_id = b)
  })
}

boot_H <- bind_rows(lapply(levels(dat_H$size_class), resample_one_size, n_boot = n_boot)) %>%
  filter(!is.na(frac_old))

sum_H <- boot_H %>%
  group_by(size_class, kingdom) %>%
  summarise(
    mean_frac = mean(frac_old),
    lower = quantile(frac_old, 0.025),
    upper = quantile(frac_old, 0.975),
    .groups = "drop"
  )

## one-sided bootstrap p: plants > animals
pvals_H <- boot_H %>%
  group_by(size_class, boot_id) %>%
  summarise(
    diff = frac_old[kingdom == "Viridiplantae"] - frac_old[kingdom == "Metazoa"],
    .groups = "drop_last"
  ) %>%
  group_by(size_class) %>%
  summarise(
    p_value = mean(diff <= 0),
    .groups = "drop"
  ) %>%
  mutate(sig = p_to_stars(p_value))

sum_H <- sum_H %>% left_join(pvals_H, by = "size_class")

SuppFig6H <- ggplot(sum_H, aes(x = size_class, y = mean_frac, fill = kingdom)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6, show.legend = TRUE) +
  geom_errorbar(
    aes(ymin = lower, ymax = upper),
    position = position_dodge(width = 0.7),
    width = 0.2,
    linewidth = 0.35,
    show.legend = FALSE
  ) +
  geom_text(
    data = sum_H %>% filter(kingdom == "Viridiplantae"),
    aes(label = sig, y = pmin(upper + 0.05, 0.98)),
    position = position_dodge(width = 0.7),
    size = 3.2,
    vjust = 0,
    show.legend = FALSE
  ) +
  scale_x_discrete(limits = c("Small", "Medium", "Large"))+
  scale_fill_manual(values = kingdom_cols) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    x = "Inversion length class",
    y = "Resampled fraction of old inversions"
  ) +
  theme_nature()+
  theme(
    axis.title.y = element_text(size = 7),
    plot.margin     = margin(1.5, 0, 1.5, 1.5, unit = "mm"))+
  theme(
    legend.position   ="right", legend.key.size = unit(3, "mm"),
    legend.title      = element_text(size = 5, lineheight = 0.9),
    legend.text       = element_text(size = 5),
    legend.margin     = margin(0, 0, 0, -0.75, unit = "mm"),
    legend.box.margin = margin(0, 0, 0, -0.75, unit = "mm")
  )


## =========================
## 4) Assemble Main Figure 3 (a–h)
## =========================
## The tag order (a–h) follows reading order left-to-right, top-to-bottom.
SuppFig6 <- ((SuppFig6A | SuppFig6B) /
           (SuppFig6C | SuppFig6D | SuppFig6E) /
           (SuppFig6F | SuppFig6G | SuppFig6H)) +
  plot_annotation(
    tag_levels = "a",
    theme = theme(plot.tag = element_text(face = "bold", size = 8))
  ) 

## Main exports (Nature-style square-ish multi-panel)
save_pdf(SuppFig6, "SuppFig6.pdf", 180, 180)
save_png(SuppFig6, "SuppFig6.png", 180, 180, dpi = 300)


## =========================
## 5) Reproducibility
## =========================
writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))