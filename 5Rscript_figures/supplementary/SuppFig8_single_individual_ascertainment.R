# single-individual ascertainment of inversion polymorphisms
#
# Input: tables.xlsx
# Required packages: readxl, dplyr, ggplot2, tidyr, purrr, patchwork

library(readxl)
library(dplyr)
library(ggplot2)
library(tidyr)
library(purrr)
library(patchwork)

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_flag <- '--file='
  path <- sub(file_flag, '', args[grepl(file_flag, args)])
  if (length(path) > 0) dirname(normalizePath(path[1])) else getwd()
}

script_dir <- get_script_dir()
PROJ_ROOT <- if (basename(script_dir) %in% c('extended_data','supplementary','main_figures')) dirname(script_dir) else script_dir
source(file.path(PROJ_ROOT, 'helpers', 'figure_plot_helpers.R'))
source(file.path(PROJ_ROOT, 'helpers', 'supp_quality_helpers.R'))

OUT_DIR <- file.path(PROJ_ROOT, 'outputs', 'supplementary')
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

infile <- "population_summary.xlsx"


dat <- read_excel(infile, sheet = "population_inversion_stats") |>
  mutate(
    kingdom = if_else(
      species %in% c(
        "Ambrosia artemisiifolia","Arabidopsis lyrata","Cenchrus americanus",
        "Helianthus annuus","Helianthus argophyllus","Helianthus petiolaris",
        "Trifolium repens","Zea mays ssp. parviglumis"
      ),
      "Plant", "Animal"
    )
  )

# IMPORTANT:
# A. lyrata is retained for population validation, but excluded from this
# detectability meta-analysis because the corrected breakpoint-marker
# reconstruction is conditioned on heterozygous marker support (HO = 1 in all
# 32 population x inversion units), so it is not an unbiased estimate of
# random-individual detectability.
d <- dat |> filter(species != "Arabidopsis lyrata")
cols <- c('Animal'='#E1812C','Plant'='#00A087')
# -------------------------
# Panel A
# -------------------------
curve <- tibble(p = seq(0, 1, length.out = 501)) |>
  mutate(HWE = 2*p*(1-p))

bins <- d |>
  mutate(bin = cut(p, breaks = seq(0,1,0.05), include.lowest = TRUE)) |>
  group_by(bin) |>
  summarise(
    p = mean(p),
    mean_HO = mean(HO),
    n = n(),
    .groups = "drop"
  )

pA <- ggplot(d, aes(p, HO)) +
  geom_point(aes(color = kingdom), alpha = 0.35, size = 1) + 
  scale_color_manual(values = cols) +
  geom_line(data = curve, aes(p, HWE), linewidth = 0.7, color = "black") +
  geom_point(data = bins, aes(p, mean_HO), size = 1.2, color = "black") +
  
  labs(
    x = "Inversion allele frequency, p",
    y = expression("Observed heterozygote frequency, " * H[O])
  ) +
  coord_cartesian(xlim = c(0,1), ylim = c(0,1)) +
  theme_classic(base_size =7) +
  theme(
    legend.position = c(0.2, 0.90), 
    legend.title = element_blank(),
    legend.text = element_text(size = 6) ,
    legend.key.size = grid::unit(3, "mm")
  )

# -------------------------
# Panel B
# -------------------------
sp <- d |>
  group_by(species, kingdom) |>
  summarise(
    mean_HO = mean(HO),
    n_units = n(),
    .groups = "drop"
  ) |>
  arrange(mean_HO) |>
  mutate(species = factor(species, levels = species))

pB <- ggplot(sp, aes(mean_HO, species, color = kingdom)) +
  geom_point(alpha = 0.65,size = 1.8)  + 
  scale_color_manual(values = cols) +
  labs(
    x = expression("Mean one-diploid detectability (" * H[O] * ")"),
    y = NULL,
    shape = NULL) +
  coord_cartesian(xlim = c(0,0.7)) +
  theme_classic(base_size =7) +
  theme(
    axis.text.y = element_text(face = "italic"), 
    legend.position = c(0.2, 0.90),  
    legend.title = element_blank(),
    legend.text = element_text(size = 6) ,
    legend.key.size = grid::unit(3, "mm")
  )



# -------------------------
# Panel C
# Equal-weight each species so Human/Helianthus do not dominate.
# For a population x inversion unit with observed HO = h,
# Pr(at least one heterozygote among n random diploids) = 1-(1-h)^n.
# -------------------------
n_grid <- 1:20

recovery <- d |>
  select(species, population, inversion_id, HO) |>
  crossing(n_diploids = n_grid) |>
  mutate(recovery = 1 - (1-HO)^n_diploids) |>
  group_by(species, n_diploids) |>
  summarise(recovery = mean(recovery), .groups = "drop") |>
  group_by(n_diploids) |>
  summarise(
    mean_recovery = mean(recovery),
    se = sd(recovery)/sqrt(n()),
    low = pmax(0, mean_recovery - 1.96*se),
    high = pmin(1, mean_recovery + 1.96*se),
    .groups = "drop"
  )

pC <- ggplot(recovery, aes(n_diploids, mean_recovery)) +
  geom_ribbon(aes(ymin = low, ymax = high), alpha = 0.2) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 1.2) +
  scale_x_continuous(breaks = c(1,2,5,10,15,20)) +
  coord_cartesian(ylim = c(0,1)) +
  labs(
    x = "Number of sampled individuals",
    y = "Probability of detecting heterozygote"
  ) +
  theme_classic(base_size = 7)


pA_ready <- pA + labs(tag = "a") + theme(plot.tag = element_text(face = 'bold', size = 8))
pB_ready <- pB + labs(tag = "b") + theme(plot.tag = element_text(face = 'bold', size = 8))
pC_ready <- pC + labs(tag = "c") + theme(plot.tag = element_text(face = 'bold', size = 8))
fig <- pA_ready | pB_ready | pC_ready


save_pdf(fig, "SuppFig8_single_individual_ascertainment.pdf", 180,60)
save_png(fig, "SuppFig8_single_individual_ascertainment.png",180,60)
