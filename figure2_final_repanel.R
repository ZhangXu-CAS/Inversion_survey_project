############################################################
## Figure 2 (Main) + Supplementary Figures S2–S3
## Inversion number & length across kingdoms
##
## Inputs (CSV):
##   - per_genome.csv      (one row per sample/genome)
##   - per_inversion.csv   (one row per inversion)
##
## Notes:
##   - Statistical comparisons: pairwise Wilcoxon rank-sum tests with BH adjustment.
##   - Bootstrap robustness: Metazoa vs Viridiplantae with equal n per kingdom.
##
## This script is a re-panel / style harmonization pass only.
## It preserves the original analytical content as much as possible.
############################################################


  library(tidyverse)
  library(ggplot2)
  library(scales)
  library(ggpubr)     # stat_pvalue_manual
  library(patchwork)  # panel assembly
  library(stringr)
  library(dplyr)
  library(readr)



## =========================
## 0) User-editable settings
## =========================
IN_PER_GENOME    <- "per_genome.csv"
IN_PER_INVERSION <- "per_inversion.csv"
OUT_DIR          <- "."

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

## Colors (consistent across panels)
kingdom_cols <- c(
  "Fungi"         = "#3274A1",
  "Metazoa"       = "#E1812C",
  "Viridiplantae" = "#00A087"
)

theme_nature <- function(base_size = 8) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid      = element_blank(),
      axis.text       = element_text(color = "black"),
      axis.title      = element_text(color = "black"),
      legend.position = "none",
      panel.border    = element_rect(color = "black", linewidth = 0.5),
      axis.ticks      = element_line(color = "black"),
      plot.title      = element_blank(),
  #    plot.subtitle   = element_blank(),
      plot.tag        = element_text(face = "bold", size = 8)
    )
}

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
p_to_star <- function(p) {
  dplyr::case_when(
    is.na(p)   ~ "",
    p < 1e-3   ~ "***",
    p < 1e-2   ~ "**",
    p < 5e-2   ~ "*",
    TRUE       ~ "ns"
  )
} 
## ---------- helpers ----------
EPS <- 1e-6
safe_log10 <- function(x, eps = EPS) {
  x <- as.numeric(x)
  log10(pmax(x, eps))
}

## Create x-axis labels with per-panel n (default)
make_x_labels_with_n <- function(df, group = "kingdom") {
  df %>%
    count(.data[[group]], name = "n") %>%
    mutate(lab = paste0(as.character(.data[[group]]), "\n(", n, ")")) %>%
    select(!!sym(group), lab) %>%
    deframe()
}

## Plain labels (no n) for panel(s) where requested
make_x_labels_plain <- function(levels_vec) {
  setNames(levels_vec, levels_vec)
}

## Pairwise Wilcoxon with BH adjustment; returns group1/group2 as strings
manual_pairwise_wilcox <- function(df, y_var, group_var = "kingdom") {
  df2 <- df %>% filter(!is.na(.data[[y_var]]), is.finite(.data[[y_var]]), !is.na(.data[[group_var]]))
  groups <- sort(unique(as.character(df2[[group_var]])))

  if (length(groups) < 2) {
    return(tibble(group1 = character(), group2 = character(),
                  p = double(), p.adj = double(), p_label = character()))
  }

  combs <- t(combn(groups, 2))

  res_list <- apply(combs, 1, function(g) {
    g1 <- g[1]; g2 <- g[2]
    y1 <- df2[df2[[group_var]] == g1, y_var, drop = TRUE]
    y2 <- df2[df2[[group_var]] == g2, y_var, drop = TRUE]
    if (sum(is.finite(y1)) < 2 || sum(is.finite(y2)) < 2) {
      tibble(group1 = g1, group2 = g2, p = NA_real_)
    } else {
      tibble(group1 = g1, group2 = g2, p = suppressWarnings(wilcox.test(y1, y2)$p.value))
    }
  })

  res <- bind_rows(res_list)

  if (all(is.na(res$p))) {
    return(res %>% mutate(p.adj = NA_real_, p_label = NA_character_))
  }

  res %>%
    mutate(
      p.adj = p.adjust(p, method = "BH"),
      p_label = case_when(
        is.na(p.adj)  ~ NA_character_,
        p.adj < 0.001 ~ "p < 0.001***",
        TRUE          ~ paste0("p = ", signif(p.adj, 2), " ", p_to_star(p.adj))
      )
    )
}

## Assign y-positions for p-value brackets
add_y_positions <- function(p_tab, max_y, step_frac = 0.08, start_frac = 0.05) {
  if (nrow(p_tab) == 0) return(p_tab)
  p_use <- p_tab %>% filter(!is.na(p_label))
  if (nrow(p_use) == 0) return(p_tab)
  if (!is.finite(max_y)) max_y <- 1

  p_use %>%
    arrange(group1, group2) %>%
    mutate(
      y.position = max_y * (1 + start_frac) + (row_number() - 1) * (max_y * step_frac)
    )
}

maybe_add_pvalues <- function(p, p_tab) {
  if (is.null(p_tab) || nrow(p_tab) == 0) return(p)
  if (!("y.position" %in% names(p_tab))) return(p)

  p_use <- p_tab %>% filter(!is.na(p_label), is.finite(y.position))
  if (nrow(p_use) == 0) return(p)

  p + stat_pvalue_manual(
    data        = p_use,
    label       = "p_label",
    y.position  = "y.position",
    xmin        = "group1",
    xmax        = "group2",
    tip.length  = 0.01,
    size        = 1.8,
    inherit.aes = FALSE
  )
}

## Unified box+jitter panel builder (kingdom-level)
boxpanel_kingdom <- function(df, y, ylab, x_labels, add_p = TRUE) {
  df2 <- df %>% filter(!is.na({{ y }}), is.finite({{ y }}))
  ycol <- rlang::as_name(rlang::enquo(y))
  ptab <- manual_pairwise_wilcox(df2, ycol) %>%
    add_y_positions(max(df2[[ycol]], na.rm = TRUE))

  p <- ggplot(df2, aes(x = kingdom, y = {{ y }}, fill = kingdom, color = kingdom)) +
    geom_boxplot(width = 0.6, alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.15, size = 1.2, alpha = 0.7) +
    scale_fill_manual(values = kingdom_cols) +
    scale_color_manual(values = kingdom_cols) +
    scale_x_discrete(labels = x_labels) +
    labs(x = NULL, y = ylab) +
    theme_nature()

  if (add_p) p <- maybe_add_pvalues(p, ptab)
  p
}

## =========================
## 1) Load and pre-process
## =========================

per_genome <- readr::read_csv(IN_PER_GENOME, show_col_types = FALSE) %>%
  filter(!is.na(kingdom)) %>%
  mutate(
    kingdom = factor(kingdom, levels = c("Fungi", "Metazoa", "Viridiplantae"))
  )

per_inv <- readr::read_csv(IN_PER_INVERSION, show_col_types = FALSE) %>%
  filter(!is.na(kingdom), !is.na(inv_len_bp)) %>%
  mutate(
    kingdom    = factor(kingdom, levels = c("Fungi", "Metazoa", "Viridiplantae")),
    inv_len_bp = as.numeric(inv_len_bp),
    log10_len  = safe_log10(inv_len_bp)
  ) %>%
  filter(is.finite(log10_len))

## S/M/L proportions (per-genome) if component counts exist
needed_sml <- c("inv_count_small_lt10kb", "inv_count_medium_10kb_1Mb", "inv_count_large_gt1Mb", "inv_count_total")
if (all(needed_sml %in% names(per_genome))) {
  per_genome <- per_genome %>%
    mutate(
      prop_small  = if_else(!is.na(inv_count_total) & inv_count_total > 0,
                            inv_count_small_lt10kb / inv_count_total, NA_real_),
      prop_medium = if_else(!is.na(inv_count_total) & inv_count_total > 0,
                            inv_count_medium_10kb_1Mb / inv_count_total, NA_real_),
      prop_large  = if_else(!is.na(inv_count_total) & inv_count_total > 0,
                            inv_count_large_gt1Mb / inv_count_total, NA_real_)
    )
}

## ---------- log10 transforms for requested metrics ----------
per_genome <- per_genome %>%
  mutate(
    inv_count_per_Gb = as.numeric(inv_count_per_Gb),
    inv_count_total  = as.numeric(inv_count_total),

    log10_inv_count_per_Gb = if_else(is.finite(inv_count_per_Gb) & inv_count_per_Gb > 0,
                                     safe_log10(inv_count_per_Gb), NA_real_),
    log10_inv_count_total  = if_else(is.finite(inv_count_total) & inv_count_total > 0,
                                     safe_log10(inv_count_total),  NA_real_),
  )



## =========================
## 2) Figure 2 (Main): 3 x 3 (a–i)
## =========================

## a) Inversion count per Gb (log10)
P2_a_dat <- per_genome %>% filter(!is.na(log10_inv_count_per_Gb))
P2_a_xlab <- make_x_labels_with_n(P2_a_dat)
P2_a <- boxpanel_kingdom(
  P2_a_dat,
  log10_inv_count_per_Gb,
  ylab = expression(paste("Inversion count per Gb (log"[10],")")),
  x_labels = P2_a_xlab
)

## b) Total inversion count (log10)
P2_b_dat <- per_genome %>% filter(!is.na(log10_inv_count_total))
P2_b_xlab <- make_x_labels_with_n(P2_b_dat)
P2_b <- boxpanel_kingdom(
  P2_b_dat,
  log10_inv_count_total,
  ylab = expression(paste("Total inversion count (log"[10],")")),
  x_labels = P2_b_xlab
)

## c) Median inversion length per genome (log10 bp)
P2_c_dat <- per_genome %>% filter(!is.na(median_log10_inv_length))
P2_c_xlab <- make_x_labels_with_n(P2_c_dat)
P2_c <- boxpanel_kingdom(
  P2_c_dat,
  median_log10_inv_length,
  ylab = expression(paste("Median inversion length (log"[10], " bp)")),
  x_labels = P2_c_xlab
)

## d) Per-inversion log10 length KDE
P2_d <- ggplot(per_inv, aes(x = log10_len, color = kingdom, fill = kingdom)) +
  geom_density(alpha = 0.25, adjust = 1.0) +
  scale_color_manual(values = kingdom_cols) +
  scale_fill_manual(values  = kingdom_cols) +
  labs(x = expression(paste("Inversion length (log"[10], " bp)")), y = "Density") +
  theme_nature() +
  theme(
    legend.position = c(0.98, 0.98),
    legend.justification = c(1, 1),
    legend.direction = "vertical",
    legend.title = element_blank(),
    legend.text = element_text(size = 6),
    legend.key.size = unit(3, "mm")
  )+
  theme(plot.margin     = margin(1.5, 3, 1.5, 2, unit = "mm"))

## e) Mean proportion of S/M/L inversions (legend inside top; no extra space; no kingdom n)
P2_e <- ggplot() + theme_void()

if (all(c("prop_small", "prop_medium", "prop_large") %in% names(per_genome))) {
  size_summary <- per_genome %>%
    filter(!is.na(prop_small), !is.na(prop_medium), !is.na(prop_large)) %>%
    pivot_longer(cols = c(prop_small, prop_medium, prop_large),
                 names_to = "size_class", values_to = "prop") %>%
    mutate(
      size_class = recode(
        size_class,
        prop_small  = "small",
        prop_medium = "medium",
        prop_large  = "large"
      ),
      size_class = factor(size_class, levels = c("small", "medium", "large"))
    ) %>%
    group_by(kingdom, size_class) %>%
    summarise(mean_prop = mean(prop, na.rm = TRUE), .groups = "drop")

  size_cols <- c(
    "small"  = "#9ecae1",
    "medium" = "#fdae6b",
    "large"  = "#a1d99b"
  )

  ## plain x labels per request (no n)
  x_plain <- make_x_labels_plain(levels(per_genome$kingdom))

  P2_e <- ggplot(size_summary, aes(x = kingdom, y = mean_prop, fill = size_class)) +
    geom_col(width = 0.7, color = "black", linewidth = 0.2) +
    scale_fill_manual(values = size_cols) +
    scale_x_discrete(labels = x_plain) +
    scale_y_continuous(
      labels = percent_format(accuracy = 1),
    ) +
    labs(x = NULL, y = "Mean proportion of inversions") +
    theme_nature() +
    theme(
      legend.position = "right",
      legend.direction = "vertical",
      legend.title = element_blank(),
      legend.text = element_text(size = 6),
      legend.key.size = unit(3, "mm"),
      plot.margin = margin(t = 2, r = 2, b = 2, l = 2)
    ) +
    theme(plot.margin     = margin(1.5, 2, 1.5, 2, unit = "mm"))
}

## g–i) Order mean summaries (log10; style matches a–c)
if (!("order" %in% names(per_genome))) {
  stop("per_genome.csv must include an 'order' column for Fig2 g–i.")
}

ord_dat <- per_genome %>%
  mutate(order = as.character(order)) %>%
  filter(!is.na(order), order != "") %>%
  group_by(kingdom, order) %>%
  summarise(
    order_mean_log10_inv_perGb  = mean(log10_inv_count_per_Gb, na.rm = TRUE),
    order_mean_median_log10_len = mean(median_log10_inv_length, na.rm = TRUE),
    order_mean_prop_large = mean(prop_large, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(kingdom = factor(kingdom, levels = c("Fungi", "Metazoa", "Viridiplantae")))


####update liner regrssion with pi --02-26-2026
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
                               method_stats = "pearson",    # "pearson" or "lm_slope"
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
# S2 panels
S2_0 <- corr_panel_kingdom(
  per_genome,
  xvar = "pi_genome_wide",
  yvar = "inv_count_per_Gb",
  xlab = "Genome-wide nucleotide diversity (π)",
  ylab = expression(paste("Inversion count per Gb (log"[10],")")),
  logy = TRUE,              
  method_stats = "pearson",
  legend_pos = c(0.70, 0.20),
  legend_text_size = 5
)

# 1) inversion counts vs π
S2_1 <- corr_panel_kingdom(
  per_genome,
  xvar = "pi_genome_wide",
  yvar = "inv_count_total",
  xlab = "Genome-wide nucleotide diversity (π)",
  ylab = expression(paste("Total inversion count (log"[10],")")),
  logy = TRUE,
  legend_pos = "right",
  legend_text_size = 5
)

# 2) median inversion length vs π（median_log10_inv_length）
S2_2 <- corr_panel_kingdom(
  per_genome,
  xvar = "pi_genome_wide",
  yvar = "median_log10_inv_length",
  xlab = "Genome-wide nucleotide diversity (π)",
  ylab = expression(paste("Median inversion length (log"[10], " bp)")),
  logy = !("median_log10_inv_length" %in% names(per_genome)),
  legend_pos = c(0.70, 0.80),
  legend_text_size = 5
)
## b) Mean inversion length per genome (log10 bp)
S2_b_dat <- per_genome %>% filter(!is.na(mean_log10_inv_length))
S2_b_xlab <- make_x_labels_with_n(S2_b_dat)
S2_b <- boxpanel_kingdom(
  S2_b_dat,
  mean_log10_inv_length,
  ylab = expression(paste("Mean inversion length (log"[10], " bp)")), 
  x_labels = S2_b_xlab
)
## Assemble Fig2 (3 x 3) in requested order

row3 <- (P2_d | P2_e | plot_spacer()) +
  plot_layout(widths = c(1.5, 1, 1e-6), guides = "keep")

Fig2 <- (( P2_a |P2_b | S2_0) /
           ( P2_c |S2_b | S2_2)/
           row3 ) +
  plot_annotation(tag_levels = "a")

save_pdf(Fig2, "Fig2.pdf", 180, 180)
save_png(Fig2, "Fig2.png", 180, 180, dpi = 300)

## =========================
## 3) Supplementary Figure S2: 2 x 2 (a–d)
## =========================

## a) Inversion count per 10,000 genes (log10)
S2_a <- ggplot() + theme_void()
if ("log10_inv_per_10k_genes" %in% names(per_genome)) {
  S2_a_dat <- per_genome %>% filter(!is.na(log10_inv_per_10k_genes))
  S2_a_xlab <- make_x_labels_with_n(S2_a_dat)
  S2_a <- boxpanel_kingdom(
    S2_a_dat,
    log10_inv_per_10k_genes,
    ylab =  expression(paste("Inversion count per 10,000 genes (log"[10],")")),
    x_labels = S2_a_xlab
  )
}

## c) Per-genome proportion of small inversions (<10 kb)
S2_c <- ggplot() + theme_void()
if ("prop_small" %in% names(per_genome)) {
  S2_c_dat <- per_genome %>% filter(!is.na(prop_small))
  S2_c_xlab <- make_x_labels_with_n(S2_c_dat)
  S2_c <- boxpanel_kingdom(
    S2_c_dat,
    prop_small,
    ylab = "Proportion of small inversions",
    x_labels = S2_c_xlab
  ) +
    scale_y_continuous(labels = percent_format(accuracy = 1))
}

## d) Per-genome proportion of medium inversions (log10)
S2_d <- ggplot() + theme_void()
if ("prop_medium" %in% names(per_genome)) {
  S2_d_dat <- per_genome %>% filter(!is.na(prop_medium))
  S2_d_xlab <- make_x_labels_with_n(S2_d_dat)
  S2_d <- boxpanel_kingdom(
    S2_d_dat,
    prop_medium,
    ylab = "Proportion of medium inversions",
    x_labels = S2_d_xlab
  )+
    scale_y_continuous(labels = percent_format(accuracy = 1))
}
## e) Per-genome proportion of large inversions (>1 Mb) (log10)
S2_e_dat <- per_genome %>% filter(!is.na(prop_large))
S2_e_xlab <- make_x_labels_with_n(S2_e_dat)
S2_e <- boxpanel_kingdom(
  S2_e_dat,
  prop_large,
  ylab = "Proportion of large inversions",
  x_labels = S2_e_xlab
)+scale_y_continuous(labels = percent_format(accuracy = 1))

## g) Order mean inversion count per Gb (log10)log10_inv_count_per_Gb
P2_g_dat <- ord_dat %>% filter(!is.na(order_mean_log10_inv_perGb), is.finite(order_mean_log10_inv_perGb))
P2_g_xlab <- make_x_labels_with_n(P2_g_dat)
P2_g <- boxpanel_kingdom(
  P2_g_dat,
  order_mean_log10_inv_perGb,
  ylab = expression(paste("Inversion count per Gb (log"[10],")")),
  x_labels = P2_g_xlab)
## h) Order mean of per-genome median inversion length (log10 bp)
P2_h_dat <- ord_dat %>% filter(!is.na(order_mean_median_log10_len), is.finite(order_mean_median_log10_len))
P2_h_xlab <- make_x_labels_with_n(P2_h_dat)
P2_h <- boxpanel_kingdom(
  P2_h_dat,
  order_mean_median_log10_len,
  ylab = expression(paste("Median inversion length (log"[10], " bp)")),
  x_labels = P2_h_xlab
)

## i) Order mean proportion of large inversions (>1 Mb) (log10)
P2_i_dat <- ord_dat %>% filter(!is.na(order_mean_prop_large), is.finite(order_mean_prop_large))
P2_i_xlab <- make_x_labels_with_n(P2_i_dat)
P2_i <- boxpanel_kingdom(
  P2_i_dat,
  order_mean_prop_large,
  ylab =  "Proportion of large inversions" ,
  x_labels = P2_i_xlab
)

####Assemble FigS2 
row1 <-(S2_a| S2_1|plot_spacer()) + plot_layout(widths = c(1.2, 1.3, 1e-6), guides = "keep")
FigS2 <-row1/ 
  (S2_c|S2_d|S2_f)/
           (P2_g | P2_h | P2_i)+
  plot_annotation(
    tag_levels = "a",
    theme = theme(plot.tag = element_text(face = "bold", size = 8))
  )


save_pdf(FigS2, "FigS2.pdf", 180, 180)
save_png(FigS2, "FigS2.png", 180, 180, dpi = 300)

## =========================
## 4) Supplementary Figure S3: 2 x 2 (a–d)
##     Per-genome distributions across major clades
## =========================

## defines “large” inversion for S3
LARGE_INV_BP <- 1e6

## Kingdom order (facet order)
KINGDOM_LEVELS <- c("Metazoa", "Viridiplantae")

find_ge1mb_col <- function(nms) {
  candidates <- nms[
    str_detect(nms, regex("1\\s*mb|1mb|>=\\s*1\\s*mb|ge\\s*1\\s*mb|1000\\s*kb|1e6|1_?000_?000", ignore_case = TRUE)) &
      str_detect(nms, regex("inv_count|count", ignore_case = TRUE))
  ]
  if (length(candidates) == 0) return(NA_character_)
  candidates[1]
}

make_clade_key <- function(df) {
  df %>%
    mutate(
      kingdom = factor(kingdom, levels = KINGDOM_LEVELS),
      clade_key = paste0(as.character(kingdom), ":", major_clade)
    )
}
clade_labels <- function(x) sub("^.*?:", "", x)

plot_panel_clade <- function(df, y, ylab, y_is_prop = FALSE) {
  y <- enquo(y)
  dat <- df %>% filter(!is.na(!!y), is.finite(!!y))

  p <- ggplot(dat, aes(x = clade_key, y = !!y)) +
    geom_boxplot(aes(fill = kingdom),
                 width = 0.55, linewidth = 0.35, outlier.shape = NA, alpha = 0.85) +
    geom_jitter(aes(color = kingdom),
                width = 0.15, size = 0.35, alpha = 0.15, show.legend = FALSE) +
    facet_wrap(~ kingdom, scales = "free_x", nrow = 1) +
    scale_fill_manual(values = kingdom_cols) +
    scale_color_manual(values = kingdom_cols) +
    scale_x_discrete(labels = clade_labels) +
    labs(x = NULL, y = ylab) +
    theme_nature() +
    theme(
      axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1, size = 5),
      axis.title.y = element_text(size = 6)
    )

  if (y_is_prop) {
    p <- p + scale_y_continuous(labels = percent_format(accuracy = 1))
  }
  p
}

pg <- readr::read_csv(IN_PER_GENOME, show_col_types = FALSE) %>%
  filter(kingdom %in% KINGDOM_LEVELS) %>%
  mutate(
    kingdom = factor(kingdom, levels = KINGDOM_LEVELS),
    major_clade = as.character(major_clade),
    inv_count_per_Gb = as.numeric(inv_count_per_Gb),
    inv_count_total  = as.numeric(inv_count_total),
    median_log10_inv_length  = as.numeric(median_log10_inv_length)
  ) %>%
  filter(!is.na(major_clade), major_clade != "")

## Derive prop_large_inv_1Mb if not present
ge1mb_col <- find_ge1mb_col(names(pg))

if (!is.na(ge1mb_col)) {
  pg <- pg %>%
    mutate(
      inv_count_ge_1Mb = as.numeric(.data[[ge1mb_col]]),
      prop_large_inv_1Mb = if_else(inv_count_total > 0, inv_count_ge_1Mb / inv_count_total, NA_real_)
    )
} else {
  ## derive from per-inversion
  pinv <- readr::read_csv(IN_PER_INVERSION, show_col_types = FALSE)
  if (!all(c("sample_id", "inv_len_bp") %in% names(pinv))) {
    stop("To derive >1Mb proportion in S3, per_inversion.csv must include: sample_id, inv_len_bp")
  }
  pinv_sum <- pinv %>%
    mutate(inv_len_bp = as.numeric(inv_len_bp)) %>%
    filter(!is.na(inv_len_bp), inv_len_bp > 0) %>%
    group_by(sample_id) %>%
    summarise(inv_count_ge_1Mb = sum(inv_len_bp > LARGE_INV_BP), .groups = "drop")

  pg <- pg %>%
    left_join(pinv_sum, by = "sample_id") %>%
    mutate(
      inv_count_ge_1Mb = replace_na(inv_count_ge_1Mb, 0),
      prop_large_inv_1Mb = if_else(inv_count_total > 0, inv_count_ge_1Mb / inv_count_total, NA_real_)
    )
}

pg <- pg %>%
  mutate(
    log10_inv_count_per_Gb = if_else(is.finite(inv_count_per_Gb) & inv_count_per_Gb > 0,
                                     safe_log10(inv_count_per_Gb), NA_real_),
    log10_inv_count_total  = if_else(is.finite(inv_count_total) & inv_count_total > 0,
                                     safe_log10(inv_count_total), NA_real_),
    median_log10_inv_length  = median_log10_inv_length
  )

pg <- make_clade_key(pg)

S3_a <- plot_panel_clade(
  pg,
  log10_inv_count_per_Gb,
  ylab = expression(paste("Inversion count per Gb (log"[10],")"))
)

S3_b <- plot_panel_clade(
  pg,
  log10_inv_count_total,
  ylab = expression(paste("Total inversion count per genome (log"[10],")"))
)

S3_c <- plot_panel_clade(
  pg,
  median_log10_inv_length,
  ylab = expression(paste("Median inversion length (log"[10], " bp)")) 
)

S3_d <- plot_panel_clade(
  pg,
  prop_large_inv_1Mb,
  ylab = "Proportion of large inversions",
  y_is_prop = TRUE
)

FigS3 <- (S3_a | S3_b) / (S3_c | S3_d) +
  plot_annotation(tag_levels = "a")

save_pdf(FigS3, "FigS3.pdf", 180, 120)
save_png(FigS3, "FigS3.png", 180, 120, dpi = 300)

## =========================
## 4) Supplementary Figure S4: 2 x 2 (a–d)
##     Bootstrapped mean
## =========================

set.seed(2025)

per_genome_MV <- per_genome %>%
  filter(kingdom %in% c("Metazoa", "Viridiplantae")) %>%
  droplevels()

# safety: ensure required columns exist
stopifnot("log10_inv_count_per_Gb" %in% names(per_genome_MV))
stopifnot("median_log10_inv_length" %in% names(per_genome_MV))

min_n_MV <- per_genome_MV %>%
  count(kingdom, name = "n") %>%
  summarise(min_n = min(n)) %>%
  pull(min_n)

# axis-label expressions (Nature-style)
lab_count_expr <- expression(paste("Inversion count per Gb (log"[10], ")"))
lab_len_expr   <- expression(paste("Median inversion length (log"[10], " bp)"))

B <- 1000
boot_res <- vector("list", B)

for (b in seq_len(B)) {
  tmp <- per_genome_MV %>%
    group_by(kingdom) %>%
    slice_sample(n = min_n_MV, replace = TRUE) %>%
    summarise(
      mean_log10_inv_count_per_Gb = mean(log10_inv_count_per_Gb, na.rm = TRUE),
      mean_median_log10_inv_len   = mean(median_log10_inv_length, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(iter = b)
  boot_res[[b]] <- tmp
}
boot_res <- bind_rows(boot_res)

# long form for CI / obs overlay
boot_long <- boot_res %>%
  pivot_longer(
    cols = c(mean_log10_inv_count_per_Gb, mean_median_log10_inv_len),
    names_to = "metric",
    values_to = "value"
  )

obs_means <- per_genome_MV %>%
  summarise(
    mean_log10_inv_count_per_Gb = mean(log10_inv_count_per_Gb, na.rm = TRUE),
    mean_median_log10_inv_len   = mean(median_log10_inv_length, na.rm = TRUE),
    .by = kingdom
  ) %>%
  pivot_longer(
    cols = c(mean_log10_inv_count_per_Gb, mean_median_log10_inv_len),
    names_to = "metric",
    values_to = "obs_mean"
  )

boot_ci <- boot_long %>%
  group_by(kingdom, metric) %>%
  summarise(
    ci_lower = quantile(value, 0.025, na.rm = TRUE),
    ci_upper = quantile(value, 0.975, na.rm = TRUE),
    .groups  = "drop"
  )

## ----------------
## S4e) bootstrap: log10 inv count per Gb
## ----------------
S4_e_dat <- boot_long %>% filter(metric == "mean_log10_inv_count_per_Gb")
S4_e_obs <- obs_means %>% filter(metric == "mean_log10_inv_count_per_Gb")
S4_e_ci  <- boot_ci   %>% filter(metric == "mean_log10_inv_count_per_Gb")

S4_e <- ggplot(S4_e_dat, aes(x = kingdom, y = value, fill = kingdom)) +
  geom_violin(trim = TRUE, alpha = 0.7) +
  geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.85, color = "black") +
  geom_point(data = S4_e_obs,
             aes(x = kingdom, y = obs_mean),
             size = 2.2, color = "black", inherit.aes = FALSE) +
  geom_errorbar(data = S4_e_ci,
                aes(x = kingdom, ymin = ci_lower, ymax = ci_upper),
                width = 0.15, linewidth = 0.6, color = "black", inherit.aes = FALSE) +
  scale_fill_manual(values = kingdom_cols) +
  labs(x = NULL, y = expression(paste("Bootstrapped mean: inversion count per Gb (log"[10], ")"))) +
  theme_nature() +
  theme(legend.position = "none", axis.title.y = element_text(size = 6))

## ----------------
## S4f) bootstrap: median inversion length (log10 bp)
## ----------------
S4_f_dat <- boot_long %>% filter(metric == "mean_median_log10_inv_len")
S4_f_obs <- obs_means %>% filter(metric == "mean_median_log10_inv_len")
S4_f_ci  <- boot_ci   %>% filter(metric == "mean_median_log10_inv_len")

S4_f <- ggplot(S4_f_dat, aes(x = kingdom, y = value, fill = kingdom)) +
  geom_violin(trim = TRUE, alpha = 0.7) +
  geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.85, color = "black") +
  geom_point(data = S4_f_obs,
             aes(x = kingdom, y = obs_mean),
             size = 2.2, color = "black", inherit.aes = FALSE) +
  geom_errorbar(data = S4_f_ci,
                aes(x = kingdom, ymin = ci_lower, ymax = ci_upper),
                width = 0.15, linewidth = 0.6, color = "black", inherit.aes = FALSE) +
  scale_fill_manual(values = kingdom_cols) +
  labs(x = NULL, y = expression(paste("Bootstrapped mean: median inversion length (log"[10], " bp)"))) +
  theme_nature() +
  theme(legend.position = "none", axis.title.y = element_text(size = 6))

## ----------------
## S4g–h) bootstrap difference distributions (Metazoa − Viridiplantae)
## ----------------
boot_wide <- boot_res %>%
  pivot_wider(
    id_cols = iter,
    names_from = kingdom,
    values_from = c(mean_log10_inv_count_per_Gb, mean_median_log10_inv_len)
  )

boot_diff <- boot_wide %>%
  transmute(
    iter,
    diff_log10_inv_count_per_Gb =
      mean_log10_inv_count_per_Gb_Metazoa - mean_log10_inv_count_per_Gb_Viridiplantae,
    diff_median_log10_inv_len =
      mean_median_log10_inv_len_Metazoa   - mean_median_log10_inv_len_Viridiplantae
  )

S4_g <- ggplot(boot_diff, aes(x = diff_log10_inv_count_per_Gb)) +
  geom_density(fill = "grey80", color = "black") +
  geom_vline(xintercept = 0, color = "red", linewidth = 0.7, linetype = "dashed") +
  labs(
    x = expression(paste(Delta, " count (Metazoa − Viridiplantae, log"[10], ")")),
    y = "Density"
  ) +
  theme_nature()

S4_h <- ggplot(boot_diff, aes(x = diff_median_log10_inv_len)) +
  geom_density(fill = "grey80", color = "black") +
  geom_vline(xintercept = 0, color = "red", linewidth = 0.7, linetype = "dashed") +
  labs(
    x = expression(paste(Delta, " length (Metazoa − Viridiplantae, log"[10], " bp)")),
    y = "Density"
  ) +
  theme_nature()

FigS4 <- (S4_e | S4_f) / (S4_g | S4_h) +
  plot_annotation(tag_levels = "a")

save_pdf(FigS4, "FigS4.pdf", 150, 120)
save_png(FigS4, "FigS4.png", 150, 120, dpi = 300)
message("Done. Outputs written to: ", normalizePath(OUT_DIR))
