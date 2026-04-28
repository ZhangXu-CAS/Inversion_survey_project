suppressPackageStartupMessages({
  library(tidyverse)
  library(ggpubr)
  library(patchwork)
  library(scales)
  library(stringr)
})

IN_PER_GENOME <- "per_genome.csv"
IN_PER_INVERSION <- "per_inversion.csv"
OUT_DIR <- "."

KINGDOM_LEVELS_ALL <- c("Fungi", "Metazoa", "Viridiplantae")
KINGDOM_LEVELS_MV  <- c("Metazoa", "Viridiplantae")
EPS <- 1e-6

kingdom_cols <- c(
  "Fungi"         = "#3274A1",
  "Metazoa"       = "#E1812C",
  "Viridiplantae" = "#00A087"
)

size_cols <- c(
  "small"  = "#9ecae1",
  "medium" = "#fdae6b",
  "large"  = "#a1d99b"
)

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

require_columns <- function(df, cols, df_name = deparse(substitute(df))) {
  missing_cols <- setdiff(cols, names(df))
  if (length(missing_cols) > 0) {
    stop(
      sprintf(
        "%s is missing required column(s): %s",
        df_name,
        paste(missing_cols, collapse = ", ")
      ),
      call. = FALSE
    )
  }
}

theme_nature <- function(base_size = 8) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid      = element_blank(),
      axis.text       = element_text(color = "black"),
      axis.title      = element_text(color = "black"),
      axis.ticks      = element_line(color = "black"),
      panel.border    = element_rect(color = "black", linewidth = 0.5),
      legend.position = "none",
      plot.title      = element_blank(),
      plot.tag        = element_text(face = "bold", size = 8)
    )
}

save_pdf <- function(plot, filename, width_mm, height_mm) {
  ggsave(
    filename = file.path(OUT_DIR, filename),
    plot = plot,
    width = width_mm,
    height = height_mm,
    units = "mm",
    useDingbats = FALSE
  )
}

save_png <- function(plot, filename, width_mm, height_mm, dpi = 300) {
  ggsave(
    filename = file.path(OUT_DIR, filename),
    plot = plot,
    width = width_mm,
    height = height_mm,
    units = "mm",
    dpi = dpi
  )
}

safe_log10 <- function(x, eps = EPS) {
  x <- as.numeric(x)
  log10(pmax(x, eps))
}

to_num <- function(x) {
  suppressWarnings(readr::parse_number(as.character(x)))
}

p_to_star <- function(p) {
  dplyr::case_when(
    is.na(p) ~ "",
    p < 1e-3 ~ "***",
    p < 1e-2 ~ "**",
    p < 5e-2 ~ "*",
    TRUE ~ "ns"
  )
}

p_label_vec <- function(p) {
  stars <- dplyr::case_when(
    is.na(p) ~ "",
    p < 1e-2 ~ "**",
    p < 5e-2 ~ "*",
    TRUE ~ ""
  )
  dplyr::case_when(
    is.na(p) ~ "p = NA",
    p < 1e-3 ~ "p < 0.001***",
    TRUE ~ paste0("p = ", signif(p, 3), stars)
  )
}

make_x_labels_with_n <- function(df, group = "kingdom") {
  df %>%
    count(.data[[group]], name = "n") %>%
    mutate(label = paste0(as.character(.data[[group]]), "\n(", n, ")")) %>%
    select(all_of(group), label) %>%
    deframe()
}

make_x_labels_plain <- function(levels_vec) {
  setNames(levels_vec, levels_vec)
}

manual_pairwise_wilcox <- function(df, y_var, group_var = "kingdom") {
  df2 <- df %>%
    filter(!is.na(.data[[y_var]]), is.finite(.data[[y_var]]), !is.na(.data[[group_var]]))

  groups <- sort(unique(as.character(df2[[group_var]])))
  if (length(groups) < 2) {
    return(tibble(group1 = character(), group2 = character(), p = double(), p.adj = double(), p_label = character()))
  }

  combs <- t(combn(groups, 2))
  res_list <- apply(combs, 1, function(g) {
    g1 <- g[1]
    g2 <- g[2]
    y1 <- df2[df2[[group_var]] == g1, y_var, drop = TRUE]
    y2 <- df2[df2[[group_var]] == g2, y_var, drop = TRUE]

    if (sum(is.finite(y1)) < 2 || sum(is.finite(y2)) < 2) {
      tibble(group1 = g1, group2 = g2, p = NA_real_)
    } else {
      tibble(group1 = g1, group2 = g2, p = suppressWarnings(wilcox.test(y1, y2, exact = FALSE)$p.value))
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
        is.na(p.adj) ~ NA_character_,
        p.adj < 0.001 ~ "p < 0.001***",
        TRUE ~ paste0("p = ", signif(p.adj, 2), " ", p_to_star(p.adj))
      )
    )
}

add_y_positions <- function(p_tab, max_y, step_frac = 0.08, start_frac = 0.05) {
  if (nrow(p_tab) == 0) return(p_tab)
  p_use <- p_tab %>% filter(!is.na(p_label))
  if (nrow(p_use) == 0) return(p_tab)
  if (!is.finite(max_y)) max_y <- 1

  p_use %>%
    arrange(group1, group2) %>%
    mutate(y.position = max_y * (1 + start_frac) + (row_number() - 1) * (max_y * step_frac))
}

maybe_add_pvalues <- function(p, p_tab) {
  if (is.null(p_tab) || nrow(p_tab) == 0) return(p)
  if (!("y.position" %in% names(p_tab))) return(p)

  p_use <- p_tab %>% filter(!is.na(p_label), is.finite(y.position))
  if (nrow(p_use) == 0) return(p)

  p + stat_pvalue_manual(
    data = p_use,
    label = "p_label",
    y.position = "y.position",
    xmin = "group1",
    xmax = "group2",
    tip.length = 0.01,
    size = 1.8,
    inherit.aes = FALSE
  )
}

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

make_legend_labels_r2p <- function(df, x, y, method = c("pearson", "lm_slope"), p_adjust = TRUE) {
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
                 if (method == "pearson") {
                   cor.test(.data[[x]], .data[[y]], method = "pearson")$p.value
                 } else {
                   summary(lm(.data[[y]] ~ .data[[x]]))$coefficients[2, 4]
                 },
                 NA_real_),
      .groups = "drop"
    )

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

corr_panel_kingdom <- function(df,
                               xvar, yvar,
                               xlab, ylab,
                               logx = FALSE, logy = FALSE,
                               method_stats = "pearson",
                               legend_pos = c(0.70, 0.20),
                               legend_text_size = 5,
                               point_alpha = 0.5,
                               point_size = 1.0,
                               line_width = 0.8) {
  d <- df %>%
    mutate(
      kingdom = factor(kingdom, levels = KINGDOM_LEVELS_ALL),
      .x = to_num(.data[[xvar]]),
      .y = to_num(.data[[yvar]])
    )

  if (logx) d <- d %>% mutate(.x = safe_log10(.x))
  if (logy) d <- d %>% mutate(.y = safe_log10(.y))

  d <- d %>% filter(!is.na(kingdom), is.finite(.x), is.finite(.y))
  legend_labels <- make_legend_labels_r2p(d, ".x", ".y", method = method_stats, p_adjust = TRUE)

  ggplot(d, aes(x = .x, y = .y, color = kingdom)) +
    geom_point(alpha = point_alpha, size = point_size) +
    geom_smooth(method = "lm", se = FALSE, linewidth = line_width) +
    scale_color_manual(
      values = kingdom_cols,
      labels = legend_labels,
      guide = guide_legend(title = NULL)
    ) +
    labs(x = xlab, y = ylab) +
    theme_nature() +
    theme(
      legend.position = legend_pos,
      legend.text = element_text(size = legend_text_size),
      legend.key.height = grid::unit(12, "pt"),
      legend.background = element_rect(fill = "transparent", color = NA)
    )
}

find_ge1mb_col <- function(nms) {
  candidates <- nms[
    str_detect(nms, regex("1\\s*mb|1mb|>=\\s*1\\s*mb|ge\\s*1\\s*mb|1000\\s*kb|1e6|1_?000_?000", ignore_case = TRUE)) &
      str_detect(nms, regex("inv_count|count", ignore_case = TRUE))
  ]
  if (length(candidates) == 0) return(NA_character_)
  candidates[1]
}

coerce_numeric_if_present <- function(df, cols) {
  cols_use <- intersect(cols, names(df))
  if (length(cols_use) == 0) return(df)
  df %>% mutate(across(all_of(cols_use), as.numeric))
}

derive_size_props <- function(per_genome) {
  needed <- c("inv_count_small_lt10kb", "inv_count_medium_10kb_1Mb", "inv_count_large_gt1Mb", "inv_count_total")
  if (!all(needed %in% names(per_genome))) return(per_genome)

  per_genome %>%
    mutate(
      prop_small  = if_else(!is.na(inv_count_total) & inv_count_total > 0,
                            inv_count_small_lt10kb / inv_count_total, NA_real_),
      prop_medium = if_else(!is.na(inv_count_total) & inv_count_total > 0,
                            inv_count_medium_10kb_1Mb / inv_count_total, NA_real_),
      prop_large  = if_else(!is.na(inv_count_total) & inv_count_total > 0,
                            inv_count_large_gt1Mb / inv_count_total, NA_real_)
    )
}

load_per_genome <- function(path = IN_PER_GENOME) {
  per_genome <- readr::read_csv(path, show_col_types = FALSE)
  require_columns(per_genome, c("kingdom"), "per_genome")

  numeric_cols <- c(
    "inv_count_per_Gb", "inv_count_total", "median_log10_inv_length", "mean_log10_inv_length",
    "pi_genome_wide", "log10_inv_per_10k_genes", "asm_scaffold_N50_Mb", "asm_total_len_bp",
    "BUSCO_completeness", "TE_total_pct", "gene_density", "indel_rate",
    "inv_count_small_lt10kb", "inv_count_medium_10kb_1Mb", "inv_count_large_gt1Mb"
  )

  per_genome <- per_genome %>%
    filter(!is.na(kingdom)) %>%
    mutate(kingdom = factor(kingdom, levels = KINGDOM_LEVELS_ALL)) %>%
    coerce_numeric_if_present(numeric_cols) %>%
    derive_size_props() %>%
    mutate(
      log10_inv_count_per_Gb = if_else(is.finite(inv_count_per_Gb) & inv_count_per_Gb > 0,
                                       safe_log10(inv_count_per_Gb), NA_real_),
      log10_inv_count_total = if_else(is.finite(inv_count_total) & inv_count_total > 0,
                                      safe_log10(inv_count_total), NA_real_)
    )

  per_genome
}

load_per_inversion <- function(path = IN_PER_INVERSION) {
  per_inv <- readr::read_csv(path, show_col_types = FALSE)
  require_columns(per_inv, c("kingdom", "inv_len_bp"), "per_inversion")

  per_inv %>%
    filter(!is.na(kingdom), !is.na(inv_len_bp)) %>%
    mutate(
      kingdom = factor(kingdom, levels = KINGDOM_LEVELS_ALL),
      inv_len_bp = as.numeric(inv_len_bp),
      log10_len = safe_log10(inv_len_bp)
    ) %>%
    filter(is.finite(log10_len))
}
