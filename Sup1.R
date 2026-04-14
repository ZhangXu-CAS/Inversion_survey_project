library(tidyverse)
library(ggpubr)   # stat_pvalue_manual, ggarrange
library(rstatix)  # add_significance

## -----------------------------
## 1) Read + rename
## -----------------------------
df_raw <- readr::read_csv("summary_stats_all.csv", show_col_types = FALSE) %>%
  dplyr::rename(
    strategy   = Strategy,
    sample     = Sample,
    kingdom    = Kingdom,
    simulated  = Simulated_INV_count,
    called_inv = Called_INV_count,
    tp         = TP,
    fp         = FP,
    fn         = FN,
    precision  = Precision,
    recall     = Recall,
    fdr        = FDR
  ) %>%
  mutate(
    kingdom = factor(kingdom, levels = c("Metazoa", "Viridiplantae", "Fungi"))
  )

## ---- strategy order (add HiFi reads-biased) ----

strategy_levels <- c("Small-biased", "Intermediate", "Large-biased", "Hifi reads-biased")


df <- df_raw %>%
  mutate(strategy = factor(strategy, levels = strategy_levels))

## -----------------------------
## 2) Collapse to one point per sample × strategy (match paired tests)
## -----------------------------
df1 <- df %>%
  group_by(sample, strategy) %>%
  summarise(
    kingdom = dplyr::first(kingdom),
    recall  = mean(recall, na.rm = TRUE),
    fdr     = mean(fdr, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  drop_na(strategy)

## -----------------------------
## 3) Paired Wilcoxon helper
## -----------------------------
pairs <- tibble(
  group1 = c("Small-biased", "Small-biased", "Small-biased",
             "Hifi reads-biased", "Hifi reads-biased",
             "Intermediate"),
  group2 = c("Hifi reads-biased", "Intermediate", "Large-biased",
             "Intermediate", "Large-biased",
             "Large-biased")
)

paired_wilcox_table <- function(dat, value_col, pairs_tbl, p_adjust = "BH") {
  res <- map2_dfr(pairs_tbl$group1, pairs_tbl$group2, function(g1, g2) {
    
    wide <- dat %>%
      filter(strategy %in% c(g1, g2)) %>%
      select(sample, strategy, {{ value_col }}) %>%
      tidyr::pivot_wider(names_from = strategy, values_from = {{ value_col }}) %>%
      tidyr::drop_na(all_of(c(g1, g2)))  # paired: both conditions present
    
    if (nrow(wide) < 2) {
      return(tibble(group1 = g1, group2 = g2, p = NA_real_, n = nrow(wide)))
    }
    
    p <- wilcox.test(wide[[g1]], wide[[g2]], paired = TRUE)$p.value
    tibble(group1 = g1, group2 = g2, p = p, n = nrow(wide))
  })
  
  res %>%
    mutate(p.adj = p.adjust(p, method = p_adjust)) %>%
    rstatix::add_significance("p.adj") %>%
    select(group1, group2, n, p, p.adj, p.adj.signif)
}

stat_recall <- paired_wilcox_table(df1, recall, pairs, p_adjust = "BH")
stat_fdr    <- paired_wilcox_table(df1, fdr,    pairs, p_adjust = "BH")

## -----------------------------
## 4) y.position for brackets (auto)
## -----------------------------
ymax_recall <- max(df1$recall, na.rm = TRUE)
ymax_fdr    <- max(df1$fdr,    na.rm = TRUE)

## 6 comparisons -> 6 y positions
stat_recall <- stat_recall %>%
  mutate(y.position = ymax_recall + seq(0.02, 0.02 * n(), by = 0.02))

stat_fdr <- stat_fdr %>%
  mutate(y.position = ymax_fdr + seq(0.01, 0.01 * n(), by = 0.01))

## -----------------------------
## 5) Colors + theme
## -----------------------------
kingdom_cols <- c(
  "Metazoa"       = "#E1812C",
  "Viridiplantae" = "#00A087",
  "Fungi"         = "grey50"
)

theme_nature_supp <- theme_classic(base_size = 8) +
  theme(
    legend.title     = element_blank(),
    legend.position  = "bottom",
    legend.direction = "horizontal",
    axis.title.x     = element_blank(),
    axis.text.x      = element_text(angle = 0, hjust = 0.5),
    plot.title       = element_text(face = "bold", size = 8)
  )

## -----------------------------
## Panel a: Recall (use df1!)
## -----------------------------
pA <- ggplot(df1, aes(x = strategy, y = recall)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, linewidth = 0.7,
               color = "black", fill = "white") +
  geom_jitter(aes(color = kingdom), width = 0.15, height = 0,
              size = 2.0, alpha = 0.7) +
  scale_color_manual(values = kingdom_cols) +
  coord_cartesian(ylim = c(min(df1$recall, na.rm = TRUE) - 0.02, ymax_recall + 0.02 * (nrow(pairs) + 1))) +
  ggpubr::stat_pvalue_manual(
    stat_recall,
    label = "p.adj.signif",   # 或改 "p.adj" 输出数值
    tip.length = 0.01,
    bracket.size = 0.5,
    size = 3
  ) +
  labs(y = "Recall (TP / (TP + FN))") +
  theme_nature_supp

## -----------------------------
## Panel b: FDR (use df1!)
## -----------------------------
pB <- ggplot(df1, aes(x = strategy, y = fdr)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, linewidth = 0.7,
               color = "black", fill = "white") +
  geom_jitter(aes(color = kingdom), width = 0.15, height = 0,
              size = 2.0, alpha = 0.7) +
  scale_color_manual(values = kingdom_cols) +
  coord_cartesian(ylim = c(0, ymax_fdr + 0.01 * (nrow(pairs) + 1))) +
  ggpubr::stat_pvalue_manual(
    stat_fdr,
    label = "p.adj.signif",
    tip.length = 0.01,
    bracket.size = 0.5,
    size = 3
  ) +
  labs(y = "FDR (FP / (TP + FP))") +
  theme_nature_supp +
  theme(legend.position = "none")  

## -----------------------------
## 6) Combine into Fig. S1 (2 panels)
## -----------------------------
figS1 <- ggpubr::ggarrange(
  pA, pB,
  ncol = 1, nrow = 2,
  labels = c("a", "b"),
  heights = c(1, 1)
)

ggsave("figS1_SyRI_simulation_strategy_focus.pdf", figS1, width = 80, height =120, units = "mm")
ggsave("figS1_SyRI_simulation_strategy_focus.png", figS1, width = 80, height =120, units = "mm", dpi = 300)

## Optional: write stats tables (handy for supplement)
write_csv(stat_recall, "figS1_stat_recall_paired_wilcox_BH.csv")
write_csv(stat_fdr,    "figS1_stat_fdr_paired_wilcox_BH.csv")