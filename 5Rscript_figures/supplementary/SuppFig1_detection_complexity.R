#!/usr/bin/env Rscript

# Supplementary Fig. 1 extension: SURVIVOR detection-complexity analysis
# Exact per_genome schema used here (no column guessing):
# sample_id, kingdom, TE_total_pct, asm_total_len_bp, pi_genome_wide

suppressPackageStartupMessages({
  library(tidyverse)
  library(glmmTMB)
  library(patchwork)
  library(scales)
})

set.seed(20260924)

# -----------------------------
# 0. Files / settings
# -----------------------------
BENCHMARK_FILE <- "summary_stats_all.csv"
PER_GENOME_FILE <- "per_genome.csv"   # change path only if the file is stored elsewhere
OUTDIR <- "detection_complexity"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

N_BOOT <- 1000L  # use 5000 for the final frozen analysis

strategy_levels <- c(
  "Small-biased",
  "Intermediate",
  "Large-biased",
  "Hifi reads-biased"
)

kingdom_levels <- c("Metazoa", "Viridiplantae", "Fungi")

kingdom_cols <- c(
  "Metazoa"       = "#E1812C",
  "Viridiplantae" = "#00A087",
  "Fungi"         = "#3274A1"
)

nature_theme <- theme_classic(base_size = 8) +
  theme(
    axis.text = element_text(colour = "black"),
    axis.title = element_text(colour = "black"),
    axis.line = element_line(linewidth = 0.45),
    axis.ticks = element_line(linewidth = 0.4),
    legend.title = element_blank(),
    legend.position = "bottom",
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 8),
    plot.margin = margin(4, 5, 4, 5)
  )

# -----------------------------
# 1. Read benchmark table
# -----------------------------
bench_raw <- readr::read_csv(BENCHMARK_FILE, show_col_types = FALSE)

required_benchmark_cols <- c(
  "Strategy", "Sample", "Kingdom",
  "Simulated_INV_count", "Called_INV_count",
  "TP", "FP", "FN", "Precision", "Recall", "FDR"
)
stopifnot(all(required_benchmark_cols %in% names(bench_raw)))

bench <- bench_raw %>%
  transmute(
    strategy   = factor(Strategy, levels = strategy_levels),
    sample_id  = as.character(Sample),
    kingdom_benchmark = as.character(Kingdom),
    simulated  = as.integer(Simulated_INV_count),
    called_inv = as.integer(Called_INV_count),
    tp         = as.integer(TP),
    fp         = as.integer(FP),
    fn         = as.integer(FN)
  ) %>%
  filter(!is.na(strategy)) %>%
  group_by(sample_id, strategy) %>%
  summarise(
    kingdom_benchmark = first(kingdom_benchmark),
    simulated = sum(simulated, na.rm = TRUE),
    called_inv = sum(called_inv, na.rm = TRUE),
    tp = sum(tp, na.rm = TRUE),
    fp = sum(fp, na.rm = TRUE),
    fn = sum(fn, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    recall = tp / (tp + fn),
    fdr = if_else(tp + fp > 0, fp / (tp + fp), NA_real_)
  )

# -----------------------------
# 2. Read the actual per_genome table using its exact headers
# -----------------------------
per_genome <- readr::read_csv(PER_GENOME_FILE, show_col_types = FALSE)

required_per_genome_cols <- c(
  "sample_id",
  "kingdom",
  "TE_total_pct",
  "asm_total_len_bp",
  "pi_genome_wide"
)
stopifnot(all(required_per_genome_cols %in% names(per_genome)))

meta <- per_genome %>%
  select(
    sample_id,
    kingdom,
    TE_total_pct,
    asm_total_len_bp,
    pi_genome_wide
  ) %>%
  mutate(
    sample_id = as.character(sample_id),
    kingdom = as.character(kingdom),
    TE_total_pct = as.numeric(TE_total_pct),
    asm_total_len_bp = as.numeric(asm_total_len_bp),
    pi_genome_wide = as.numeric(pi_genome_wide)
  )

if (anyDuplicated(meta$sample_id)) {
  stop("per_genome contains duplicated sample_id values; resolve before analysis.")
}


x <- bench %>%
  left_join(meta, by = "sample_id")


# Descriptive dataset includes fungi; primary inference is Metazoa vs Viridiplantae.
x <- x %>%
  mutate(
    kingdom = factor(kingdom, levels = kingdom_levels),
    log10_genome_size = log10(asm_total_len_bp),
    log10_pi = log10(pi_genome_wide)
  )

model_df <- x %>%
  filter(kingdom %in% c("Metazoa", "Viridiplantae")) %>%
  droplevels()

# Scale continuous predictors once across unique benchmark genomes.
scaling_df <- model_df %>%
  distinct(sample_id, TE_total_pct, log10_genome_size, log10_pi)

te_mu <- mean(scaling_df$TE_total_pct)
te_sd <- sd(scaling_df$TE_total_pct)
gs_mu <- mean(scaling_df$log10_genome_size)
gs_sd <- sd(scaling_df$log10_genome_size)
pi_mu <- mean(scaling_df$log10_pi)
pi_sd <- sd(scaling_df$log10_pi)

model_df <- model_df %>%
  mutate(
    TE_z = (TE_total_pct - te_mu) / te_sd,
    genome_size_z = (log10_genome_size - gs_mu) / gs_sd,
    pi_z = (log10_pi - pi_mu) / pi_sd,
    kingdom = factor(kingdom, levels = c("Metazoa", "Viridiplantae")),
    strategy = factor(strategy, levels = strategy_levels),
    sample_id = factor(sample_id)
  )

write_csv(model_df, file.path(OUTDIR, "model_input_exact.csv"))

# Correlations among complexity covariates at the genome level.
cor_dat <- model_df %>%
  distinct(sample_id, kingdom, TE_total_pct, log10_genome_size, log10_pi)
cor_mat <- cor(
  cor_dat %>% select(TE_total_pct, log10_genome_size, log10_pi),
  method = "spearman",
  use = "complete.obs"
)
write_csv(
  as.data.frame(cor_mat) %>% rownames_to_column("variable"),
  file.path(OUTDIR, "complexity_spearman_correlations.csv")
)

# -----------------------------
# 3. Main beta-binomial mixed models
# -----------------------------
fit_recall <- glmmTMB(
  cbind(tp, fn) ~ strategy + kingdom + TE_z + genome_size_z + pi_z + (1 | sample_id),
  family = betabinomial(link = "logit"),
  data = model_df
)

fit_fdr <- glmmTMB(
  cbind(fp, tp) ~ strategy + kingdom + TE_z + genome_size_z + pi_z + (1 | sample_id),
  family = betabinomial(link = "logit"),
  data = model_df
)

if (!isTRUE(fit_recall$sdr$pdHess)) warning("Recall model Hessian is not positive definite.")
if (!isTRUE(fit_fdr$sdr$pdHess)) warning("FDR model Hessian is not positive definite.")

# Save fixed-effect coefficients with Wald CIs (secondary; probability-scale effects below are primary).
extract_fixed <- function(fit, outcome) {
  cc <- summary(fit)$coefficients$cond
  tibble(
    outcome = outcome,
    term = rownames(cc),
    estimate_logit = cc[, "Estimate"],
    SE = cc[, "Std. Error"],
    z = cc[, "z value"],
    p = cc[, "Pr(>|z|)"],
    conf_low_logit = estimate_logit - 1.96 * SE,
    conf_high_logit = estimate_logit + 1.96 * SE
  )
}

main_coef <- bind_rows(
  extract_fixed(fit_recall, "Recall"),
  extract_fixed(fit_fdr, "FDR")
)
write_csv(main_coef, file.path(OUTDIR, "main_model_coefficients.csv"))

# -----------------------------
# 4. Reviewer-targeted interaction sensitivity models
#    Fit one interaction at a time to avoid over-parameterizing 20 plant/animal genomes.
# -----------------------------
interaction_models <- list(
  recall_TE_by_kingdom = glmmTMB(
    cbind(tp, fn) ~ strategy + kingdom * TE_z + genome_size_z + pi_z + (1 | sample_id),
    family = betabinomial(link = "logit"), data = model_df
  ),
  recall_size_by_kingdom = glmmTMB(
    cbind(tp, fn) ~ strategy + kingdom * genome_size_z + TE_z + pi_z + (1 | sample_id),
    family = betabinomial(link = "logit"), data = model_df
  ),
  recall_TE_by_regime = glmmTMB(
    cbind(tp, fn) ~ strategy * TE_z + kingdom + genome_size_z + pi_z + (1 | sample_id),
    family = betabinomial(link = "logit"), data = model_df
  ),
  fdr_TE_by_kingdom = glmmTMB(
    cbind(fp, tp) ~ strategy + kingdom * TE_z + genome_size_z + pi_z + (1 | sample_id),
    family = betabinomial(link = "logit"), data = model_df
  ),
  fdr_size_by_kingdom = glmmTMB(
    cbind(fp, tp) ~ strategy + kingdom * genome_size_z + TE_z + pi_z + (1 | sample_id),
    family = betabinomial(link = "logit"), data = model_df
  ),
  fdr_TE_by_regime = glmmTMB(
    cbind(fp, tp) ~ strategy * TE_z + kingdom + genome_size_z + pi_z + (1 | sample_id),
    family = betabinomial(link = "logit"), data = model_df
  )
)

interaction_coef <- imap_dfr(interaction_models, ~ extract_fixed(.x, .y))
write_csv(interaction_coef, file.path(OUTDIR, "interaction_model_coefficients.csv"))

# -----------------------------
# 5. Adjusted effects on probability scale
# -----------------------------
q_te <- quantile(scaling_df$TE_total_pct, c(0.25, 0.75), names = FALSE)
q_gs <- quantile(scaling_df$log10_genome_size, c(0.25, 0.75), names = FALSE)
q_pi <- quantile(scaling_df$log10_pi, c(0.25, 0.75), names = FALSE)

qvals <- list(
  TE_z = (q_te - te_mu) / te_sd,
  genome_size_z = (q_gs - gs_mu) / gs_sd,
  pi_z = (q_pi - pi_mu) / pi_sd
)

predict_fixed <- function(fit, nd) {
  as.numeric(predict(fit, newdata = nd, type = "response", re.form = NA, allow.new.levels = TRUE))
}

average_contrast <- function(fit, dat, variable, low = NULL, high = NULL,
                             level1 = NULL, level2 = NULL) {
  d1 <- dat
  d2 <- dat
  if (!is.null(low) && !is.null(high)) {
    d1[[variable]] <- low
    d2[[variable]] <- high
  } else {
    d1[[variable]] <- factor(level1, levels = levels(dat[[variable]]))
    d2[[variable]] <- factor(level2, levels = levels(dat[[variable]]))
  }
  mean(predict_fixed(fit, d2) - predict_fixed(fit, d1), na.rm = TRUE)
}

calculate_effects <- function(fit, ref_dat) {
  tibble(
    effect = c(
      "TE content: Q75 - Q25",
      "Genome size: Q75 - Q25",
      "Haplotypic divergence: Q75 - Q25",
      "Viridiplantae - Metazoa",
      "Intermediate - Small-biased",
      "Large-biased - Small-biased",
      "HiFi reads-biased - Small-biased"
    ),
    estimate = c(
      average_contrast(fit, ref_dat, "TE_z", qvals$TE_z[1], qvals$TE_z[2]),
      average_contrast(fit, ref_dat, "genome_size_z", qvals$genome_size_z[1], qvals$genome_size_z[2]),
      average_contrast(fit, ref_dat, "pi_z", qvals$pi_z[1], qvals$pi_z[2]),
      average_contrast(fit, ref_dat, "kingdom", level1 = "Metazoa", level2 = "Viridiplantae"),
      average_contrast(fit, ref_dat, "strategy", level1 = "Small-biased", level2 = "Intermediate"),
      average_contrast(fit, ref_dat, "strategy", level1 = "Small-biased", level2 = "Large-biased"),
      average_contrast(fit, ref_dat, "strategy", level1 = "Small-biased", level2 = "Hifi reads-biased")
    )
  )
}

point_recall <- calculate_effects(fit_recall, model_df) %>% mutate(outcome = "Recall")
point_fdr <- calculate_effects(fit_fdr, model_df) %>% mutate(outcome = "FDR")

# -----------------------------
# 6. Genome-cluster bootstrap CIs
# -----------------------------
boot_one <- function(b, dat) {
  ids <- levels(dat$sample_id)
  draw <- sample(ids, size = length(ids), replace = TRUE)
  
  boot_dat <- map2_dfr(draw, seq_along(draw), function(id, k) {
    dat %>%
      filter(as.character(sample_id) == id) %>%
      mutate(sample_boot = factor(paste0(id, "__", k)))
  })
  
  # sample_boot is the bootstrap cluster; all 4 regimes of a sampled genome stay together.
  fr <- try(
    glmmTMB(
      cbind(tp, fn) ~ strategy + kingdom + TE_z + genome_size_z + pi_z + (1 | sample_boot),
      family = betabinomial(link = "logit"), data = boot_dat
    ), silent = TRUE
  )
  
  ff <- try(
    glmmTMB(
      cbind(fp, tp) ~ strategy + kingdom + TE_z + genome_size_z + pi_z + (1 | sample_boot),
      family = betabinomial(link = "logit"), data = boot_dat
    ), silent = TRUE
  )
  
  out <- list()
  
  # Prediction reference must contain the bootstrap grouping variable even with re.form = NA.
  ref_boot <- model_df %>%
    mutate(sample_boot = factor(paste0(as.character(sample_id), "__REF")))
  
  if (!inherits(fr, "try-error") && isTRUE(fr$sdr$pdHess)) {
    out[["Recall"]] <- calculate_effects(fr, ref_boot) %>%
      mutate(outcome = "Recall", bootstrap = b)
  }
  
  if (!inherits(ff, "try-error") && isTRUE(ff$sdr$pdHess)) {
    out[["FDR"]] <- calculate_effects(ff, ref_boot) %>%
      mutate(outcome = "FDR", bootstrap = b)
  }
  
  bind_rows(out)
}

boot_res <- map_dfr(seq_len(N_BOOT), boot_one, dat = model_df)
write_csv(boot_res, file.path(OUTDIR, "cluster_bootstrap_effects.csv"))

boot_summary <- boot_res %>%
  group_by(outcome, effect) %>%
  summarise(
    n_success = n(),
    conf_low = quantile(estimate, 0.025, na.rm = TRUE),
    conf_high = quantile(estimate, 0.975, na.rm = TRUE),
    .groups = "drop"
  )

point_effects <- bind_rows(point_recall, point_fdr) %>%
  left_join(boot_summary, by = c("outcome", "effect")) %>%
  mutate(
    estimate_pp = 100 * estimate,
    conf_low_pp = 100 * conf_low,
    conf_high_pp = 100 * conf_high
  )

write_csv(point_effects, file.path(OUTDIR, "adjusted_effects_probability_scale.csv"))

# -----------------------------
# 7. Raw descriptive panels (same visual logic as original Fig. S1)
# -----------------------------
pA <- ggplot(x, aes(strategy, recall)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, linewidth = 0.6,
               colour = "black", fill = "white") +
  geom_jitter(aes(colour = kingdom_benchmark), width = 0.12, size = 1.7, alpha = 0.75) +
  scale_colour_manual(values = kingdom_cols, drop = FALSE) +
  labs(x = NULL, y = "Recall") +
  nature_theme +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))

pB <- ggplot(x, aes(strategy, fdr)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, linewidth = 0.6,
               colour = "black", fill = "white") +
  geom_jitter(aes(colour = kingdom_benchmark), width = 0.12, size = 1.7, alpha = 0.75) +
  scale_colour_manual(values = kingdom_cols, drop = FALSE) +
  labs(x = NULL, y = "False discovery rate") +
  nature_theme +
  theme(axis.text.x = element_text(angle = 25, hjust = 1), legend.position = "none")

# -----------------------------
# 8. Effect-size forest plots
# -----------------------------
effect_levels <- rev(c(
  "TE content: Q75 - Q25",
  "Genome size: Q75 - Q25",
  "Haplotypic divergence: Q75 - Q25",
  "Viridiplantae - Metazoa",
  "Intermediate - Small-biased",
  "Large-biased - Small-biased",
  "HiFi reads-biased - Small-biased"
))

plot_effects <- function(dat, outcome_name, xlab) {
  dd <- dat %>%
    filter(outcome == outcome_name) %>%
    mutate(effect = factor(effect, levels = effect_levels))
  
  ggplot(dd, aes(x = estimate_pp, y = effect)) +
    geom_vline(xintercept = 0, linetype = 2, linewidth = 0.45, colour = "grey40") +
    geom_errorbarh(aes(xmin = conf_low_pp, xmax = conf_high_pp), height = 0.16, linewidth = 0.6) +
    geom_point(size = 2.0) +
    labs(x = xlab, y = NULL) +
    nature_theme +
    theme(legend.position = "none")
}

pC <- plot_effects(
  point_effects, "Recall",
  "Adjusted change in recall (percentage points)"
)

pD <- plot_effects(
  point_effects, "FDR",
  "Adjusted change in FDR (percentage points)"
)

figS1 <- ((pA / pB) | (pC / pD)) + plot_annotation(tag_levels = "a")

ggsave(
  file.path(OUTDIR, "figS1_SURVIVOR_detection_complexity.pdf"),
  figS1, width = 175, height = 150, units = "mm"
)
ggsave(
  file.path(OUTDIR, "figS1_SURVIVOR_detection_complexity.png"),
  figS1, width = 175, height = 150, units = "mm", dpi = 600
)

# -----------------------------
# 9. Diagnostic raw relationships (not necessarily for publication)
# -----------------------------
genome_perf <- x %>%
  group_by(sample_id, kingdom_benchmark, TE_total_pct, asm_total_len_bp, pi_genome_wide) %>%
  summarise(
    mean_recall = mean(recall, na.rm = TRUE),
    mean_fdr = mean(fdr, na.rm = TRUE),
    .groups = "drop"
  )

p_te <- ggplot(genome_perf, aes(TE_total_pct, mean_recall, colour = kingdom_benchmark)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.6) +
  scale_colour_manual(values = kingdom_cols) +
  labs(x = "Genome-wide TE content (%)", y = "Mean recall") + nature_theme

p_size <- ggplot(genome_perf, aes(log10(asm_total_len_bp), mean_recall, colour = kingdom_benchmark)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.6) +
  scale_colour_manual(values = kingdom_cols) +
  labs(x = expression(log[10]~"assembly size (bp)"), y = "Mean recall") + nature_theme

p_pi <- ggplot(genome_perf, aes(log10(pi_genome_wide), mean_recall, colour = kingdom_benchmark)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.6) +
  scale_colour_manual(values = kingdom_cols) +
  labs(x = expression(log[10]~pi[genome]), y = "Mean recall") + nature_theme

p_fdr_te <- ggplot(genome_perf, aes(TE_total_pct, mean_fdr, colour = kingdom_benchmark)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.6) +
  scale_colour_manual(values = kingdom_cols) +
  labs(x = "Genome-wide TE content (%)", y = "Mean FDR") + nature_theme

diag_fig <- (p_te | p_size) / (p_pi | p_fdr_te) + plot_annotation(tag_levels = "a")

ggsave(
  file.path(OUTDIR, "diagnostic_complexity_relationships.png"),
  diag_fig, width = 170, height = 140, units = "mm",dpi = 600
)
