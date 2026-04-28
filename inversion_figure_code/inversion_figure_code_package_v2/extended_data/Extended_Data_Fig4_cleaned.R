############################################################
## Extended Data Figure 4
## Bootstrap robustness (Metazoa vs Viridiplantae)
##
## Panel order:
##   a  Bootstrapped mean inversion count per Gb
##   b  Bootstrapped mean median inversion length
##   c  Distribution of bootstrap differences in inversion count per Gb
##   d  Distribution of bootstrap differences in median inversion length
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

set.seed(2025)
B <- 1000

per_genome <- load_per_genome() %>%
  filter(kingdom %in% KINGDOM_LEVELS_MV) %>%
  droplevels()

require_columns(per_genome, c("log10_inv_count_per_Gb", "median_log10_inv_length"), "per_genome")

min_n <- per_genome %>%
  count(kingdom, name = "n") %>%
  summarise(min_n = min(n)) %>%
  pull(min_n)

boot_res <- vector("list", B)
for (b in seq_len(B)) {
  boot_res[[b]] <- per_genome %>%
    group_by(kingdom) %>%
    slice_sample(n = min_n, replace = TRUE) %>%
    summarise(
      mean_log10_inv_count_per_Gb = mean(log10_inv_count_per_Gb, na.rm = TRUE),
      mean_median_log10_inv_len = mean(median_log10_inv_length, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(iter = b)
}
boot_res <- bind_rows(boot_res)

boot_long <- boot_res %>%
  pivot_longer(
    cols = c(mean_log10_inv_count_per_Gb, mean_median_log10_inv_len),
    names_to = "metric",
    values_to = "value"
  )

obs_means <- per_genome %>%
  summarise(
    mean_log10_inv_count_per_Gb = mean(log10_inv_count_per_Gb, na.rm = TRUE),
    mean_median_log10_inv_len = mean(median_log10_inv_length, na.rm = TRUE),
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
    .groups = "drop"
  )

## a
ed4a_dat <- boot_long %>% filter(metric == "mean_log10_inv_count_per_Gb")
ed4a_obs <- obs_means %>% filter(metric == "mean_log10_inv_count_per_Gb")
ed4a_ci  <- boot_ci %>% filter(metric == "mean_log10_inv_count_per_Gb")

ed4a <- ggplot(ed4a_dat, aes(x = kingdom, y = value, fill = kingdom)) +
  geom_violin(trim = TRUE, alpha = 0.7) +
  geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.85, color = "black") +
  geom_point(data = ed4a_obs, aes(x = kingdom, y = obs_mean), size = 2.2, color = "black", inherit.aes = FALSE) +
  geom_errorbar(
    data = ed4a_ci,
    aes(x = kingdom, ymin = ci_lower, ymax = ci_upper),
    width = 0.15,
    linewidth = 0.6,
    color = "black",
    inherit.aes = FALSE
  ) +
  scale_fill_manual(values = kingdom_cols) +
  labs(x = NULL, y = expression(paste("Bootstrapped mean: inversion count per Gb (log"[10], ")"))) +
  theme_nature() +
  theme(legend.position = "none", axis.title.y = element_text(size = 6))

## b
ed4b_dat <- boot_long %>% filter(metric == "mean_median_log10_inv_len")
ed4b_obs <- obs_means %>% filter(metric == "mean_median_log10_inv_len")
ed4b_ci  <- boot_ci %>% filter(metric == "mean_median_log10_inv_len")

ed4b <- ggplot(ed4b_dat, aes(x = kingdom, y = value, fill = kingdom)) +
  geom_violin(trim = TRUE, alpha = 0.7) +
  geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.85, color = "black") +
  geom_point(data = ed4b_obs, aes(x = kingdom, y = obs_mean), size = 2.2, color = "black", inherit.aes = FALSE) +
  geom_errorbar(
    data = ed4b_ci,
    aes(x = kingdom, ymin = ci_lower, ymax = ci_upper),
    width = 0.15,
    linewidth = 0.6,
    color = "black",
    inherit.aes = FALSE
  ) +
  scale_fill_manual(values = kingdom_cols) +
  labs(x = NULL, y = expression(paste("Bootstrapped mean: median inversion length (log"[10], " bp)"))) +
  theme_nature() +
  theme(legend.position = "none", axis.title.y = element_text(size = 6))

boot_wide <- boot_res %>%
  pivot_wider(
    id_cols = iter,
    names_from = kingdom,
    values_from = c(mean_log10_inv_count_per_Gb, mean_median_log10_inv_len)
  )

boot_diff <- boot_wide %>%
  transmute(
    iter,
    diff_log10_inv_count_per_Gb = mean_log10_inv_count_per_Gb_Metazoa - mean_log10_inv_count_per_Gb_Viridiplantae,
    diff_median_log10_inv_len = mean_median_log10_inv_len_Metazoa - mean_median_log10_inv_len_Viridiplantae
  )

## c
ed4c <- ggplot(boot_diff, aes(x = diff_log10_inv_count_per_Gb)) +
  geom_density(fill = "grey80", color = "black") +
  geom_vline(xintercept = 0, color = "red", linewidth = 0.7, linetype = "dashed") +
  labs(
    x = expression(paste(Delta, " count (Metazoa - Viridiplantae, log"[10], ")")),
    y = "Density"
  ) +
  theme_nature()

## d
ed4d <- ggplot(boot_diff, aes(x = diff_median_log10_inv_len)) +
  geom_density(fill = "grey80", color = "black") +
  geom_vline(xintercept = 0, color = "red", linewidth = 0.7, linetype = "dashed") +
  labs(
    x = expression(paste(Delta, " length (Metazoa - Viridiplantae, log"[10], " bp)")),
    y = "Density"
  ) +
  theme_nature()

Extended_Data_Fig4 <- (ed4a | ed4b) / (ed4c | ed4d) +
  plot_annotation(tag_levels = "a")

save_pdf(Extended_Data_Fig4, "Extended_Data_Fig4.pdf", 150, 120)
save_png(Extended_Data_Fig4, "Extended_Data_Fig4.png", 150, 120, dpi = 300)

message("Done: Extended_Data_Fig4.pdf / Extended_Data_Fig4.png")
