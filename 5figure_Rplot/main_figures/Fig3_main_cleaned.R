
get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_flag <- '--file='
  path <- sub(file_flag, '', args[grepl(file_flag, args)])
  if (length(path) > 0) dirname(normalizePath(path[1])) else getwd()
}
SCRIPT_DIR <- get_script_dir()
PROJ_ROOT <- normalizePath(file.path(SCRIPT_DIR, '..'), winslash = '/', mustWork = FALSE)
source(file.path(PROJ_ROOT, 'helpers', 'figure_plot_helpers.R'))
source(file.path(PROJ_ROOT, 'helpers', 'figure34_helpers.R'))
OUT_DIR <- file.path(PROJ_ROOT, 'outputs', 'main_figures')
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
IN_PER_GENOME <- file.path(PROJ_ROOT, 'per_genome.csv')
IN_PER_INVERSION <- file.path(PROJ_ROOT, 'per_inversion.csv')


############################################################
## Figure 3 (main)
############################################################
AGE_THRESHOLD <- 1.0
MAX_POINTS_PER_KINGDOM <- 10000L
SEED <- 123
per_inv <- load_per_inversion_age(IN_PER_INVERSION, age_threshold = AGE_THRESHOLD)
per_inv_age <- per_inv %>% filter(!is.na(pi_ratio), is.finite(log10_len))
log10_min <- max(2, floor(min(per_inv_age$log10_len, na.rm = TRUE)))
log10_max <- min(8, ceiling(max(per_inv_age$log10_len, na.rm = TRUE)))
age_max_plot <- quantile(per_inv_age$pi_ratio, 0.99, na.rm = TRUE)
set.seed(SEED)
plot_dat <- per_inv_age
smooth_dat <- per_inv_age %>% filter(log10_len >= 3, log10_len <= 8) %>% group_by(kingdom) %>% group_modify(~{n_here <- nrow(.x); n_take <- min(MAX_POINTS_PER_KINGDOM, n_here); if (n_here > n_take) slice_sample(.x, n = n_take) else .x}) %>% ungroup()
fig3a <- ggplot(plot_dat, aes(x = log10_len, y = pi_ratio)) + geom_point(aes(color = kingdom), alpha = 0.10, size = 0.3, show.legend = FALSE) + geom_smooth(data = smooth_dat, aes(color = kingdom), method = 'loess', se = TRUE, linewidth = 0.8, span = 0.7, show.legend = FALSE) + scale_color_manual(values = kingdom_cols) + facet_wrap(~ kingdom, nrow = 1) + coord_cartesian(ylim = c(0, age_max_plot)) + geom_hline(yintercept = AGE_THRESHOLD, linetype = 'dashed', linewidth = 0.35) + labs(x = expression(paste('Inversion length (log'[10], ' bp)')), y = expression(pi[inversion] / pi[genome])) + theme_nature()
n_len_bins <- 50; n_age_bins <- 50
len_breaks <- seq(log10_min, log10_max, length.out = n_len_bins + 1)
age_breaks <- seq(0, age_max_plot, length.out = n_age_bins + 1)
dat_grid <- per_inv_age %>% mutate(len_bin = cut(log10_len, breaks = len_breaks, include.lowest = TRUE, labels = FALSE), age_bin = cut(pi_ratio, breaks = age_breaks, include.lowest = TRUE, labels = FALSE)) %>% drop_na(len_bin, age_bin) %>% count(kingdom, len_bin, age_bin, name = 'n') %>% group_by(kingdom) %>% mutate(prop = n/sum(n), prop_rel = prop/max(prop), prop_plot = sqrt(prop_rel), len_mid = (len_breaks[len_bin] + len_breaks[len_bin+1])/2, age_mid = (age_breaks[age_bin] + age_breaks[age_bin+1])/2) %>% ungroup()
fig3b <- ggplot(dat_grid, aes(x = len_mid, y = age_mid)) + geom_tile(aes(fill = prop_plot)) + geom_hline(yintercept = AGE_THRESHOLD, linetype = 'dashed', linewidth = 0.35, colour = 'black') + facet_wrap(~ kingdom, nrow = 1) + scale_fill_gradientn(colours = c('#f7fbff','#deebf7','#9ecae1','#3182bd','#31a354','#ffff00'), values = rescale(c(0,0.12,0.28,0.5,0.75,1)), limits = c(0,1), name = 'Relative density', guide = guide_colorbar(barheight = grid::unit(15,'mm'), barwidth = grid::unit(3,'mm'), ticks.colour = 'black', frame.colour = 'black')) + labs(x = expression(paste('Inversion length (log'[10], ' bp)')), y = expression(pi[inversion] / pi[genome])) + theme_nature() + theme(legend.title = element_text(size = 5, lineheight = 0.9), legend.text = element_text(size = 5), legend.key.size = grid::unit(3,'mm'))
glm_dat <- per_inv_age %>% filter(!is.na(age_class), log10_len >= log10_min, log10_len <= log10_max)
glm_old <- glm(I(age_class == 'Old') ~ log10_len * kingdom, data = glm_dat, family = binomial())
newdat <- expand.grid(log10_len = seq(log10_min, log10_max, length.out = 200), kingdom = c('Metazoa','Viridiplantae'))
pred <- predict(glm_old, newdata = newdat, se.fit = TRUE)
newdat <- newdat %>% mutate(fit = pred$fit, se = pred$se.fit, p = plogis(fit), p_low = plogis(fit - 1.96 * se), p_high = plogis(fit + 1.96 * se))
fig3c <- ggplot(newdat, aes(x = log10_len, y = p, colour = kingdom, fill = kingdom)) + geom_ribbon(aes(ymin = p_low, ymax = p_high), alpha = 0.2, colour = NA, show.legend = TRUE) + geom_line(linewidth = 0.9, show.legend = TRUE) + scale_colour_manual(values = kingdom_cols) + scale_fill_manual(values = kingdom_cols) + coord_cartesian(ylim = c(0,1)) + labs(x = expression(paste('Inversion length (log'[10], ' bp)')), y = 'Predicted P(old)') + theme_nature() + theme(legend.position = c(0.75,0.2), legend.background = element_rect(fill = NA), legend.key.size = grid::unit(3,'mm'), legend.text = element_text(size = 5))
bin_width <- 0.5; n_bins <- ceiling((log10_max - log10_min)/bin_width)
dat_d <- per_inv_age %>% filter(log10_len >= log10_min, log10_len <= log10_max, !is.na(age_class)) %>% mutate(bin_index = pmin(pmax(floor((log10_len - log10_min)/bin_width), 0L), n_bins - 1L), len_mid = log10_min + bin_width * (bin_index + 0.5)) %>% group_by(kingdom, len_mid) %>% summarise(n_total = n(), n_old = sum(age_class == 'Old'), frac_old = n_old/n_total, .groups = 'drop')
fig3d <- ggplot(dat_d, aes(x = len_mid, y = frac_old, color = kingdom)) + geom_line(linewidth = 0.9, show.legend = TRUE) + geom_point(size = 1.4, show.legend = TRUE) + scale_color_manual(values = kingdom_cols) + coord_cartesian(ylim = c(0,1)) + labs(x = expression(paste('Inversion length (log'[10], ' bp, binned)')), y = 'Fraction of old inversions') + theme_nature() + theme(legend.position = c(0.75,0.2), legend.background = element_rect(fill = NA), legend.key.size = grid::unit(3,'mm'), legend.text = element_text(size = 5))
dat_e <- per_inv_age %>% filter(!is.na(size_class), log10_len >= log10_min, log10_len <= log10_max)
pi99_e <- quantile(dat_e$pi_ratio, 0.99, na.rm = TRUE)
stat_e <- dat_e %>% group_by(size_class) %>% summarise(p_value = tryCatch(wilcox.test(pi_ratio ~ kingdom, exact = FALSE)$p.value, error = function(e) NA_real_), .groups = 'drop') %>% mutate(sig = p_to_star(p_value))
pos_e <- dat_e %>% group_by(size_class) %>% summarise(y = quantile(pi_ratio, 0.98, na.rm = TRUE), .groups = 'drop') %>% mutate(y = pmin(y + 0.06*pi99_e, pi99_e*0.98))
stat_e <- left_join(stat_e, pos_e, by = 'size_class')
fig3e <- ggplot(dat_e, aes(x = size_class, y = pi_ratio, fill = kingdom)) + geom_point(aes(color = kingdom), position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.7), alpha = 0.25, size = 0.45, show.legend = TRUE) + geom_boxplot(position = position_dodge(width = 0.7), width = 0.45, outlier.shape = NA, linewidth = 0.35, show.legend = TRUE) + geom_hline(yintercept = AGE_THRESHOLD, linetype = 'dashed', linewidth = 0.35) + geom_text(data = stat_e, aes(x = size_class, y = y, label = sig), inherit.aes = FALSE, size = 3.2, vjust = 0) + scale_color_manual(values = kingdom_cols) + scale_fill_manual(values = kingdom_cols) + coord_cartesian(ylim = c(0, pi99_e)) + labs(x = 'Inversion length class', y = expression(pi[inversion] / pi[genome])) + theme_nature() + theme(plot.margin = margin(1.5,0,1.5,1.5, unit='mm'), legend.position='right', legend.key.size = grid::unit(3,'mm'), legend.text = element_text(size = 5))
dat_f <- per_inv_age %>% filter(inv_len_bp >= 1e4, log10_len >= log10_min, log10_len <= log10_max)
pi95_f <- quantile(dat_f$pi_ratio, 0.95, na.rm = TRUE)
fig3f <- ggplot(dat_f, aes(x = pi_ratio, colour = kingdom)) + stat_ecdf(aes(y = after_stat(1 - y)), linewidth = 0.7, show.legend = TRUE) + scale_colour_manual(values = kingdom_cols) + coord_cartesian(xlim = c(0, pi95_f)) + labs(x = expression(pi[inversion] / pi[genome]), y = expression(paste('Fraction of inversions >= ', pi[inversion]/pi[genome]))) + theme_nature() + theme(legend.position = c(0.75,0.8), axis.title.y = element_text(size = 7), legend.background = element_rect(fill = NA), legend.key.size = grid::unit(3,'mm'), legend.text = element_text(size = 5))
threshold_grid <- seq(0.3,1.5,by=0.1)
dat_large <- per_inv_age %>% filter(size_class == 'Large', log10_len >= log10_min, log10_len <= log10_max)
sens_res <- purrr::map_dfr(threshold_grid, function(thr) dat_large %>% mutate(is_old_thr = pi_ratio > thr) %>% group_by(kingdom) %>% summarise(n_total = n(), frac_old = if_else(n_total > 0, mean(is_old_thr, na.rm = TRUE), NA_real_), .groups = 'drop') %>% mutate(age_threshold = thr))
fig3g <- ggplot(sens_res, aes(x = age_threshold, y = frac_old, colour = kingdom)) + geom_line(linewidth = 0.9, show.legend = TRUE) + geom_point(size = 1.2, show.legend = TRUE) + scale_colour_manual(values = kingdom_cols) + coord_cartesian(ylim = c(0,1)) + labs(x = 'Different thresholds for old', y = 'Fraction of old inversions (>=1 Mb)') + theme_nature() + theme(legend.position = c(0.3,0.2), axis.title.y = element_text(size = 7), legend.background = element_rect(fill = NA), legend.key.size = grid::unit(3,'mm'), legend.text = element_text(size = 5))
dat_h <- per_inv_age %>% filter(!is.na(size_class), !is.na(age_class))
resample_one_size <- function(size_label, n_boot = 1000) { dat_sc <- dat_h %>% filter(size_class == size_label); tab_sc <- table(dat_sc$kingdom); if (length(tab_sc) < 2) return(NULL); n_min <- min(tab_sc); purrr::map_dfr(seq_len(n_boot), function(b) { dat_b <- dat_sc %>% group_by(kingdom) %>% slice_sample(n = n_min, replace = FALSE) %>% ungroup(); dat_b %>% group_by(kingdom) %>% summarise(frac_old = mean(age_class == 'Old'), .groups = 'drop') %>% mutate(size_class = size_label, boot_id = b) }) }
boot_h <- bind_rows(lapply(levels(dat_h$size_class), resample_one_size, n_boot = 1000)) %>% filter(!is.na(frac_old))
sum_h <- boot_h %>% group_by(size_class, kingdom) %>% summarise(mean_frac = mean(frac_old), lower = quantile(frac_old, 0.025), upper = quantile(frac_old, 0.975), .groups = 'drop')
pvals_h <- boot_h %>% group_by(size_class, boot_id) %>% summarise(diff = frac_old[kingdom == 'Viridiplantae'] - frac_old[kingdom == 'Metazoa'], .groups = 'drop_last') %>% group_by(size_class) %>% summarise(p_value = mean(diff <= 0), .groups = 'drop') %>% mutate(sig = p_to_star(p_value))
sum_h <- left_join(sum_h, pvals_h, by = 'size_class')
fig3h <- ggplot(sum_h, aes(x = size_class, y = mean_frac, fill = kingdom)) + geom_col(position = position_dodge(width = 0.7), width = 0.6, show.legend = TRUE) + geom_errorbar(aes(ymin = lower, ymax = upper), position = position_dodge(width = 0.7), width = 0.2, linewidth = 0.35, show.legend = FALSE) + geom_text(data = sum_h %>% filter(kingdom == 'Viridiplantae'), aes(label = sig, y = pmin(upper + 0.05,0.98)), position = position_dodge(width = 0.7), size = 3.2, vjust = 0, show.legend = FALSE) + scale_x_discrete(limits = c('Small','Medium','Large')) + scale_fill_manual(values = kingdom_cols) + coord_cartesian(ylim = c(0,1)) + labs(x = 'Inversion length class', y = 'Resampled fraction of old inversions') + theme_nature() + theme(axis.title.y = element_text(size = 7), plot.margin = margin(1.5,0,1.5,1.5, unit='mm'), legend.position='right', legend.key.size = grid::unit(3,'mm'), legend.text = element_text(size = 5))
Fig3 <- ((fig3a | fig3b) / (fig3c | fig3d | fig3e) / (fig3f | fig3g | fig3h)) + plot_annotation(tag_levels = 'a', theme = theme(plot.tag = element_text(face = 'bold', size = 8)))
save_pdf(Fig3, 'Fig3.pdf', 180, 180); save_png(Fig3, 'Fig3.png', 180, 180, dpi = 300); writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, 'Fig3_sessionInfo.txt'))
