
get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_flag <- '--file='
  path <- sub(file_flag, '', args[grepl(file_flag, args)])
  if (length(path) > 0) dirname(normalizePath(path[1])) else getwd()
}
SCRIPT_DIR <- get_script_dir()
PROJ_ROOT <- normalizePath(file.path(SCRIPT_DIR, '..'), winslash = '/', mustWork = FALSE)
source(file.path(PROJ_ROOT, 'helpers', 'figure_plot_helpers.R'))
source(file.path(PROJ_ROOT, 'helpers', 'supp_quality_helpers.R'))
OUT_DIR <- file.path(PROJ_ROOT, 'outputs', 'supplementary')
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
IN_PER_GENOME <- file.path(PROJ_ROOT, 'per_genome.csv')
IN_PER_INVERSION <- file.path(PROJ_ROOT, 'per_inversion.csv')


############################################################
## Supplementary Figure 1
############################################################
IN_SUMMARY <- file.path(PROJ_ROOT, 'summary_stats_all.csv')
if (!file.exists(IN_SUMMARY)) stop('Supplementary Fig. 1 expects summary_stats_all.csv in the package root.', call. = FALSE)
df_raw <- readr::read_csv(IN_SUMMARY, show_col_types = FALSE) %>% dplyr::rename(strategy = Strategy, sample = Sample, kingdom = Kingdom, simulated = Simulated_INV_count, called_inv = Called_INV_count, tp = TP, fp = FP, fn = FN, precision = Precision, recall = Recall, fdr = FDR) %>% mutate(kingdom = factor(kingdom, levels = c('Metazoa','Viridiplantae','Fungi')))
strategy_levels <- c('Small-biased','Intermediate','Large-biased','Hifi reads-biased')
df <- df_raw %>% mutate(strategy = factor(strategy, levels = strategy_levels))
df1 <- df %>% group_by(sample, strategy) %>% summarise(kingdom = dplyr::first(kingdom), recall = mean(recall, na.rm = TRUE), fdr = mean(fdr, na.rm = TRUE), .groups = 'drop') %>% drop_na(strategy)
pairs <- tibble(group1 = c('Small-biased','Small-biased','Small-biased','Hifi reads-biased','Hifi reads-biased','Intermediate'), group2 = c('Hifi reads-biased','Intermediate','Large-biased','Intermediate','Large-biased','Large-biased'))
paired_wilcox_table <- function(dat, value_col, pairs_tbl, p_adjust = 'BH') { res <- purrr::map2_dfr(pairs_tbl$group1, pairs_tbl$group2, function(g1,g2) { wide <- dat %>% filter(strategy %in% c(g1,g2)) %>% select(sample, strategy, {{value_col}}) %>% tidyr::pivot_wider(names_from = strategy, values_from = {{value_col}}) %>% tidyr::drop_na(all_of(c(g1,g2))); if (nrow(wide) < 2) return(tibble(group1 = g1, group2 = g2, p = NA_real_, n = nrow(wide))); tibble(group1 = g1, group2 = g2, p = wilcox.test(wide[[g1]], wide[[g2]], paired = TRUE, exact = FALSE)$p.value, n = nrow(wide)) }); res %>% mutate(p.adj = p.adjust(p, method = p_adjust)) %>% rstatix::add_significance('p.adj') %>% select(group1,group2,n,p,p.adj,p.adj.signif) }
stat_recall <- paired_wilcox_table(df1, recall, pairs); stat_fdr <- paired_wilcox_table(df1, fdr, pairs)
ymax_recall <- max(df1$recall, na.rm = TRUE); ymax_fdr <- max(df1$fdr, na.rm = TRUE)
stat_recall <- stat_recall %>% mutate(y.position = ymax_recall + seq(0.02, 0.02*n(), by = 0.02)); stat_fdr <- stat_fdr %>% mutate(y.position = ymax_fdr + seq(0.01, 0.01*n(), by = 0.01))
theme_nature_supp <- theme_classic(base_size = 8) + theme(legend.title = element_blank(), legend.position = 'bottom', legend.direction = 'horizontal', axis.title.x = element_blank(), axis.text.x = element_text(angle = 0, hjust = 0.5))
sim_cols <- c('Metazoa'='#E1812C','Viridiplantae'='#00A087','Fungi'='grey50')
pA <- ggplot(df1, aes(x = strategy, y = recall)) + geom_boxplot(width = 0.55, outlier.shape = NA, linewidth = 0.7, color = 'black', fill = 'white') + geom_jitter(aes(color = kingdom), width = 0.15, height = 0, size = 2.0, alpha = 0.7) + scale_color_manual(values = sim_cols) + coord_cartesian(ylim = c(min(df1$recall, na.rm = TRUE)-0.02, ymax_recall + 0.02*(nrow(pairs)+1))) + ggpubr::stat_pvalue_manual(stat_recall, label = 'p.adj.signif', tip.length = 0.01, bracket.size = 0.5, size = 3) + labs(y = 'Recall (TP / (TP + FN))') + theme_nature_supp
pB <- ggplot(df1, aes(x = strategy, y = fdr)) + geom_boxplot(width = 0.55, outlier.shape = NA, linewidth = 0.7, color = 'black', fill = 'white') + geom_jitter(aes(color = kingdom), width = 0.15, height = 0, size = 2.0, alpha = 0.7) + scale_color_manual(values = sim_cols) + coord_cartesian(ylim = c(0, ymax_fdr + 0.01*(nrow(pairs)+1))) + ggpubr::stat_pvalue_manual(stat_fdr, label = 'p.adj.signif', tip.length = 0.01, bracket.size = 0.5, size = 3) + labs(y = 'FDR (FP / (TP + FP))') + theme_nature_supp + theme(legend.position = 'none')
SuppFig1 <- ggpubr::ggarrange(pA, pB, ncol = 1, nrow = 2, labels = c('a','b'), heights = c(1,1))
save_pdf(SuppFig1, 'SuppFig1_SyRI_simulation_strategy.pdf', 80, 120); save_png(SuppFig1, 'SuppFig1_SyRI_simulation_strategy.png', 80, 120, dpi = 300); readr::write_csv(stat_recall, file.path(OUT_DIR, 'SuppFig1_stat_recall_paired_wilcox_BH.csv')); readr::write_csv(stat_fdr, file.path(OUT_DIR, 'SuppFig1_stat_fdr_paired_wilcox_BH.csv')); writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, 'SuppFig1_sessionInfo.txt'))
