############################################################
## Supplementary Figures: Robustness to genome-quality filters
##
## Goal:
##   Recompute/visualize Fig. 2A–C + all Fig. S2 panels under
##   genome-quality filters for reproducibility/robustness.
##
## Filters:
##   Supp. Fig. 2: BUSCO_completeness >= 80
##   Supp. Fig. 3: BUSCO_completeness >= 80 AND BUSCO_duplication <= 50
##   Supp. Fig. 4: LOCO robustness (random 50 samples Insects)
## Inputs:
##   - per_genome.csv
##   - per_inversion.csv
##
############################################################

  library(tidyverse)
  library(ggplot2)
  library(scales)
  library(ggpubr)
  library(patchwork)

## =========================
## 0) User-editable settings
## =========================
IN_PER_GENOME    <- "per_genome.csv"
IN_PER_INVERSION <- "per_inversion.csv"
OUT_DIR          <- "suppfig_busco_outputs"

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

## Colors
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
      axis.text.x     = element_text(color = "black", lineheight = 0.9, margin = margin(t = 2, unit = "mm")),
      axis.title      = element_text(color = "black"),
      legend.position = "none",
      panel.border    = element_rect(color = "black", linewidth = 0.5),
      axis.ticks      = element_line(color = "black"),
      plot.title      = element_blank(),
      plot.subtitle   = element_blank(),
     axis.title.y = element_text(size = 7)
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
make_x_labels_n <- function(df) {
  ## Two-line labels to prevent crowding: Kingdom \n (n)
  df %>%
    count(kingdom, name = "n") %>%
    mutate(lab = paste0(as.character(kingdom), "\n(", n, ")")) %>%
    select(kingdom, lab) %>%
    deframe()
}
EPS <- 1e-6
safe_log10 <- function(x, eps = EPS) {
  x <- as.numeric(x)
  log10(pmax(x, eps))
}

manual_pairwise_wilcox <- function(df, y_var, group_var = "kingdom") {
  df2 <- df %>% filter(!is.na(.data[[y_var]]), !is.na(.data[[group_var]]))
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
    if (sum(!is.na(y1)) < 2 || sum(!is.na(y2)) < 2) {
      tibble(group1 = g1, group2 = g2, p = NA_real_)
    } else {
      tibble(group1 = g1, group2 = g2, p = wilcox.test(y1, y2)$p.value)
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
        p.adj < 0.001 ~ "p < 0.001",
        TRUE          ~ paste0("p = ", signif(p.adj, 2))
      )
    )
}

add_y_positions <- function(p_tab, max_y) {
  if (nrow(p_tab) == 0) return(p_tab)
  p_use <- p_tab %>% filter(!is.na(p_label))
  if (nrow(p_use) == 0) return(p_tab)
  if (!is.finite(max_y)) max_y <- 1

  p_use %>%
    arrange(group1, group2) %>%
    mutate(y.position = seq(from = max_y * 1.05, to = max_y * 1.25, length.out = n()))
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

## =========================
## 1) Load
## =========================
per_genome_raw <- readr::read_csv(IN_PER_GENOME, show_col_types = FALSE) %>%
  filter(!is.na(kingdom)) %>%
  mutate(
    kingdom = factor(kingdom, levels = c("Fungi", "Metazoa", "Viridiplantae")),
    BUSCO_completeness = readr::parse_number(as.character(BUSCO_completeness)),
    BUSCO_duplication  = readr::parse_number(as.character(BUSCO_duplication))
  )

per_inv_raw <- readr::read_csv(IN_PER_INVERSION, show_col_types = FALSE) %>%
  filter(!is.na(kingdom), !is.na(inv_len_bp)) %>%
  mutate(
    kingdom    = factor(kingdom, levels = c("Fungi", "Metazoa", "Viridiplantae")),
    inv_len_bp = as.numeric(inv_len_bp),
    log10_len  = log10(inv_len_bp)
  ) %>%
  filter(is.finite(log10_len))

## Helper: compute S/M/L proportions if component counts exist
add_sml_props <- function(pg) {
  needed_sml <- c("inv_count_small_lt10kb", "inv_count_medium_10kb_1Mb", "inv_count_large_gt1Mb", "inv_count_total")
  if (all(needed_sml %in% names(pg))) {
    pg %>%
      mutate(
        prop_small  = if_else(!is.na(inv_count_total) & inv_count_total > 0,
                              inv_count_small_lt10kb / inv_count_total, NA_real_),
        prop_medium = if_else(!is.na(inv_count_total) & inv_count_total > 0,
                              inv_count_medium_10kb_1Mb / inv_count_total, NA_real_),
        prop_large  = if_else(!is.na(inv_count_total) & inv_count_total > 0,
                              inv_count_large_gt1Mb / inv_count_total, NA_real_)
      )
  } else {
    pg
  }
}

## =========================
## 2) Filtered datasets
## =========================
make_filtered <- function(filter_name = c("BUSCO80", "BUSCO80_DUP50")) {
  filter_name <- match.arg(filter_name)

  pg <- per_genome_raw %>%
    filter(!is.na(BUSCO_completeness), BUSCO_completeness >= 80)

  if (filter_name == "BUSCO80_DUP50") {
    pg <- pg %>%
      filter(!is.na(BUSCO_duplication), BUSCO_duplication <= 50)
  }

  pg <- add_sml_props(pg)

  keep_samples <- unique(pg$sample_id)

  pi <- per_inv_raw %>%
    filter(sample_id %in% keep_samples)

  list(per_genome = pg, per_inv = pi)
}

ds1 <- make_filtered("BUSCO80")
ds2 <- make_filtered("BUSCO80_DUP50")

## Save filtered tables for transparency
readr::write_csv(ds1$per_genome, file.path(OUT_DIR, "per_genome.BUSCO80.csv"))
readr::write_csv(ds2$per_genome, file.path(OUT_DIR, "per_genome.BUSCO80_DUP50.csv"))
readr::write_csv(ds1$per_inv,    file.path(OUT_DIR, "per_inversion.BUSCO80.csv"))
readr::write_csv(ds2$per_inv,    file.path(OUT_DIR, "per_inversion.BUSCO80_DUP50.csv"))

## =========================
## 3) Panel builders (reusable)
## =========================
panel_box <- function(df, y, ylab, filename = NULL, width_mm = 70, height_mm = 65) {
  df2 <- df %>% filter(!is.na(.data[[y]]))
  if (nrow(df2) == 0) return(ggplot() + theme_void() + labs(title = "No data"))

  xlab <- make_x_labels_n(df2)

  ptab <- manual_pairwise_wilcox(df2, y) %>%
    add_y_positions(max(df2[[y]], na.rm = TRUE))

  p <- ggplot(df2, aes(x = kingdom, y = .data[[y]], fill = kingdom, color = kingdom)) +
    geom_boxplot(width = 0.6, alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.15, size = 1.2, alpha = 0.7) +
    scale_fill_manual(values = kingdom_cols) +
    scale_color_manual(values = kingdom_cols) +
    scale_x_discrete(labels = xlab) +
    labs(x = NULL, y = ylab) +
    theme_nature()

  p <- maybe_add_pvalues(p, ptab)

  if (!is.null(filename)) save_pdf(p, filename, width_mm, height_mm)
  p
}

panel_density_len <- function(pi, filename = NULL, width_mm = 70, height_mm = 65, add_legend = TRUE) {
  if (nrow(pi) == 0) return(ggplot() + theme_void() + labs(title = "No inversions"))

  log10_min <- floor(quantile(pi$log10_len, 0.01, na.rm = TRUE))
  log10_max <- ceiling(quantile(pi$log10_len, 0.99, na.rm = TRUE))

  p <- ggplot(pi, aes(x = log10_len, color = kingdom, fill = kingdom)) +
    geom_density(alpha = 0.25, adjust = 1.0) +
    scale_color_manual(values = kingdom_cols) +
    scale_fill_manual(values  = kingdom_cols) +
    coord_cartesian(xlim = c(log10_min, log10_max)) +
    labs(x = expression(paste("Inversion length (log"[10], " bp)")), y = "Density") +
    theme_nature()

  if (isTRUE(add_legend)) {
    p <- p + theme(
      legend.position = c(0.98, 0.98),
      legend.justification = c(1, 1),
      legend.direction = "vertical",
      legend.title = element_blank(),
      legend.text = element_text(size = 6),
      legend.key.size = unit(3, "mm"),
      legend.background = element_rect(fill = alpha("white", 0.75), color = "black", linewidth = 0.2),
      legend.margin = margin(1, 1, 1, 1, unit = "mm")
    )
  }

  if (!is.null(filename)) save_pdf(p, filename, width_mm, height_mm)
  p
}

panel_inv_per_10k_genes <- function(pg, filename = NULL, width_mm = 70, height_mm = 65) {
  if (!("gene_counts" %in% names(pg))) {
    return(ggplot() + theme_void() + labs(title = "S2B unavailable\n(gene_counts missing)"))
  }

  df <- pg %>%
    filter(!is.na(inv_count_total), !is.na(gene_counts)) %>%
    mutate(
      inv_count_total_num = readr::parse_number(as.character(inv_count_total)),
      gene_counts_num     = readr::parse_number(as.character(gene_counts))
    ) %>%
    filter(!is.na(inv_count_total_num), !is.na(gene_counts_num), gene_counts_num > 0) %>%
    mutate(inv_per_10k_genes_log10 = log10(inv_count_total_num / gene_counts_num * 10000))

  if (nrow(df) == 0) {
    return(ggplot() + theme_void() + labs(title = "S2B unavailable\n(no valid rows)"))
  }

  xlab <- make_x_labels_n(df)
  ptab <- manual_pairwise_wilcox(df, "inv_per_10k_genes_log10") %>%
    add_y_positions(max(df$inv_per_10k_genes_log10, na.rm = TRUE))

  p <- ggplot(df, aes(x = kingdom, y = inv_per_10k_genes_log10, fill = kingdom, color = kingdom)) +
    geom_boxplot(width = 0.6, alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.15, size = 1.2, alpha = 0.7) +
    scale_fill_manual(values = kingdom_cols) +
    scale_color_manual(values = kingdom_cols) +
    scale_x_discrete(labels = xlab) +
    labs(x = NULL, y =  expression(paste("Inversion count per 10,000 genes (log"[10],")"))) +
    theme_nature()

  p <- maybe_add_pvalues(p, ptab)

  if (!is.null(filename)) save_pdf(p, filename, width_mm, height_mm)
  p
}

panel_prop <- function(pg, prop_col, ylab, filename = NULL, width_mm = 70, height_mm = 65) {
  if (!(prop_col %in% names(pg))) {
    return(ggplot() + theme_void() + labs(title = paste0(ylab, "\nunavailable")))
  }

  df <- pg %>% filter(!is.na(.data[[prop_col]]))
  if (nrow(df) == 0) {
    return(ggplot() + theme_void() + labs(title = paste0(ylab, "\nunavailable")))
  }

  xlab <- make_x_labels_n(df)
  ptab <- manual_pairwise_wilcox(df, prop_col) %>%
    add_y_positions(max(df[[prop_col]], na.rm = TRUE))

  p <- ggplot(df, aes(x = kingdom, y = .data[[prop_col]], fill = kingdom, color = kingdom)) +
    geom_hline(yintercept = 0, color = "grey85", linewidth = 0.3) +
    geom_boxplot(width = 0.6, alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.15, size = 1.2, alpha = 0.7) +
    scale_fill_manual(values = kingdom_cols) +
    scale_color_manual(values = kingdom_cols) +
    scale_x_discrete(labels = xlab) +
    coord_cartesian(ylim = c(0, 1.3)) +
    labs(x = NULL, y = ylab) +
    theme_nature()

  p <- maybe_add_pvalues(p, ptab)

  if (!is.null(filename)) save_pdf(p, filename, width_mm, height_mm)
  p
}

## =========================
## 4) Build two 9-panel supplementary figures
## =========================
build_supp9 <- function(ds, prefix) {
  pg <- ds$per_genome %>%
    mutate(
      inv_count_per_Gb = as.numeric(inv_count_per_Gb),
      inv_count_total  = as.numeric(inv_count_total),
      
      log10_inv_count_per_Gb = if_else(is.finite(inv_count_per_Gb) & inv_count_per_Gb > 0,
                                       safe_log10(inv_count_per_Gb), NA_real_),
      log10_inv_count_total  = if_else(is.finite(inv_count_total) & inv_count_total > 0,
                                       safe_log10(inv_count_total),  NA_real_),
    )
  
  pi <- ds$per_inv

  ## Row 1
  P1 <- panel_box(pg, "log10_inv_count_per_Gb", expression(paste("Inversion count per Gb (log"[10],")")), filename = paste0(prefix, "_A.pdf"))
  P2 <- panel_box(pg, "log10_inv_count_total",  expression(paste("Total inversion count (log"[10],")")), filename = paste0(prefix, "_B.pdf"))
  P3 <- panel_inv_per_10k_genes(pg, filename = paste0(prefix, "_C.pdf"))

  ## Row 2
  P4 <- panel_box(pg, "median_log10_inv_length",
                  expression(paste("Median inversion length (log"[10], " bp)")),
                  filename = paste0(prefix, "_D.pdf"))
  P5 <- panel_box(pg, "mean_log10_inv_length",
                  expression(paste("Mean inversion length (log"[10], " bp)")),
                  filename = paste0(prefix, "_E.pdf"))
  P6 <- panel_density_len(pi, filename = paste0(prefix, "_F.pdf"), add_legend = TRUE)

  ## Row 3
  P7 <- panel_prop(pg, "prop_small",  "Proportion of small inversions",        filename = paste0(prefix, "_G.pdf"))
  P8 <- panel_prop(pg, "prop_medium", "Proportion of medium inversions",  filename = paste0(prefix, "_H.pdf"))
  P9 <- panel_prop(pg, "prop_large",  "Proportion of large inversions",        filename = paste0(prefix, "_I.pdf"))

  Fig <- (P1 | P2 | P3) / (P4 | P5 | P6) / (P7 | P8 | P9) +
    plot_annotation(
      tag_levels = "a",
      theme = theme(plot.tag = element_text(face = "bold", size = 8))
    )

  save_pdf(Fig, paste0(prefix, ".pdf"), 180, 170)
  save_png(Fig, paste0(prefix, ".png"), 180, 170)
  invisible(Fig)
}

## Supplementary Figure 2: BUSCO >= 80
build_supp9(ds1, "SuppFig2_BUSCO80")

## Supplementary Figure 3: BUSCO >= 80 AND duplication <= 50
build_supp9(ds2, "SuppFig3_BUSCO80_dup50")

## Session info
writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))

## =========================
##  Supp Fig 4
## =========================

OUT_DIR_LOCO <- "suppfig_LOCO_noInsects_outputs"
dir.create(OUT_DIR_LOCO, showWarnings = FALSE, recursive = TRUE)

save_pdf_loco <- function(plot, filename, width_mm, height_mm) {
  ggsave(
    filename = file.path(OUT_DIR_LOCO, filename),
    plot     = plot,
    width    = width_mm,
    height   = height_mm,
    units    = "mm",
    useDingbats = FALSE
  )
}
save_png_loco <- function(plot, filename, width_mm, height_mm) {
  ggsave(
    filename = file.path(OUT_DIR_LOCO, filename),
    plot     = plot,
    width    = width_mm,
    height   = height_mm,
    units    = "mm",
    dpi      = 300
  )
}


build_supp9_loco <- function(ds, prefix) {
  pg <- ds$per_genome %>%
    mutate(
      inv_count_per_Gb = as.numeric(inv_count_per_Gb),
      inv_count_total  = as.numeric(inv_count_total),
      
      log10_inv_count_per_Gb = if_else(is.finite(inv_count_per_Gb) & inv_count_per_Gb > 0,
                                       safe_log10(inv_count_per_Gb), NA_real_),
      log10_inv_count_total  = if_else(is.finite(inv_count_total) & inv_count_total > 0,
                                       safe_log10(inv_count_total),  NA_real_),
    )
  
  pi <- ds$per_inv
  
  ## Row 1
  P1 <- panel_box(pg, "log10_inv_count_per_Gb", expression(paste("Inversion count per Gb (log"[10],")")), filename = paste0(prefix, "_A.pdf"))
  P2 <- panel_box(pg, "log10_inv_count_total",  expression(paste("Total inversion count (log"[10],")")), filename = paste0(prefix, "_B.pdf"))
  P3 <- panel_inv_per_10k_genes(pg, filename = paste0(prefix, "_C.pdf"))
  
  ## Row 2
  P4 <- panel_box(pg, "median_log10_inv_length",
                  expression(paste("Median inversion length (log"[10], " bp)")),
                  filename = paste0(prefix, "_D.pdf"))
  P5 <- panel_box(pg, "mean_log10_inv_length",
                  expression(paste("Mean inversion length (log"[10], " bp)")),
                  filename = paste0(prefix, "_E.pdf"))
  P6 <- panel_density_len(pi, filename = paste0(prefix, "_F.pdf"), add_legend = TRUE)
  
  ## Row 3
  P7 <- panel_prop(pg, "prop_small",  "Proportion of small inversions",
                   filename = paste0(prefix, "_G.pdf"))
  P8 <- panel_prop(pg, "prop_medium", "Proportion of medium inversions",
                   filename = paste0(prefix, "_H.pdf"))
  P9 <- panel_prop(pg, "prop_large",  "Proportion of large inversions",
                   filename = paste0(prefix, "_I.pdf"))
  
  Fig <- (P1 | P2 | P3) / (P4 | P5 | P6) / (P7 | P8 | P9) +
    plot_annotation(
      tag_levels = "a",
      theme = theme(plot.tag = element_text(face = "bold", size = 8))
    )
  
  save_pdf_loco(Fig, paste0(prefix, ".pdf"), 180, 170)
  save_png_loco(Fig, paste0(prefix, ".png"), 180, 170)
  invisible(Fig)
}

set.seed(20260127)

ds_loco <- {
  pg0 <- per_genome_raw %>%
    mutate(major_clade = as.character(major_clade))
  
  insects_ids <- pg0 %>%
    filter(kingdom == "Metazoa", !is.na(major_clade), major_clade == "Insects") %>%
    distinct(sample_id) %>%
    pull(sample_id)
  
  n_keep <- min(50L, length(insects_ids))
  
  insects_keep <- if (n_keep > 0) sample(insects_ids, size = n_keep, replace = FALSE) else character(0)
  
  ## insects_keep
  pg <- pg0 %>%
    filter(
      !(kingdom == "Metazoa" & !is.na(major_clade) & major_clade == "Insects") |
        sample_id %in% insects_keep
    ) %>%
    add_sml_props()
  
  keep_samples <- unique(pg$sample_id)
  
  pi <- per_inv_raw %>%
    filter(sample_id %in% keep_samples)
  
  list(per_genome = pg, per_inv = pi)
}


readr::write_csv(ds_loco$per_genome, file.path(OUT_DIR_LOCO, "per_genome.LOCO_noInsects.csv"))
readr::write_csv(ds_loco$per_inv,    file.path(OUT_DIR_LOCO, "per_inversion.LOCO_noInsects.csv"))

## (5)9 panel（= Supplementary Fig. 4）
build_supp9_loco(ds_loco, "SuppFig4_LOCO_noInsects")

## sessionInfo
writeLines(capture.output(sessionInfo()), file.path(OUT_DIR_LOCO, "sessionInfo.txt"))
