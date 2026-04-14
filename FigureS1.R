## =========================
## Figure 1: Global overview of >1,000 haplotype-resolved genomes
## =========================

library(tidyverse)
library(ggridges)
library(patchwork)
library(ggpubr)
library(scales)
library(rstatix)
library(ggplot2)
library(dplyr)

## -------------------------
## Load data
## -------------------------

pg <- read_csv("per_genome.csv", show_col_types = FALSE) %>%
  filter(kingdom %in% c("Metazoa","Viridiplantae","Fungi")) %>%
  mutate(
    kingdom = factor(kingdom, levels=c("Fungi","Metazoa","Viridiplantae")),
    TE_total_pct = as.numeric(TE_total_pct),
    gene_density = as.numeric(gene_density),
    pi_genome_wide = as.numeric(pi_genome_wide),
    indel_rate = as.numeric(indel_rate),
    asm_scaffold_N50_Mb = as.numeric(asm_scaffold_N50_Mb),
    asm_total_len_bp = as.numeric(asm_total_len_bp)
  )

## -------------------------
## Color palette & theme
## -------------------------

kingdom_cols <- c(
  "Fungi"="#3274A1",
  "Metazoa"="#E1812C",
  "Viridiplantae"="#00A087"
)

theme_nature <- function(base_size=8){
  theme_bw(base_size=base_size) +
    theme(
      panel.grid = element_blank(),
      axis.text = element_text(color="black"),
      axis.title = element_text(color="black"),
      panel.border = element_rect(color="black", linewidth=0.5),
      legend.position="none"
    )
}
## -------------------------
## pairwise comparisons & helper theme
## -------------------------
cmp_pairs <- list(
  c("Fungi","Metazoa"),
  c("Fungi","Viridiplantae"),
  c("Metazoa","Viridiplantae")
)

p_label <- function(p){
  stars <- dplyr::case_when(
    p < 1e-3 ~ "***",  
    p < 1e-2 ~ "**",
    p < 5e-2 ~ "*",
    TRUE     ~ ""
  )
  
  dplyr::if_else(
    p < 1e-3,
    "p < 0.001***",
    paste0("p = ", signif(p, 3), stars)
  )
}

pg <- pg %>%
  mutate(kingdom = factor(kingdom, levels = c("Fungi","Metazoa","Viridiplantae")))
## =========================
## 1A: sampling overview
## =========================

p1A <- ggplot(pg, aes(x=kingdom, fill=kingdom))+
  geom_bar(color="black", width=0.7)+
  scale_fill_manual(values=kingdom_cols)+
  labs(y="Number of genomes", x=NULL, title="")+
  theme_nature()
## =========================
## 1B: Assembly size vs N50
## =========================

p1B <- ggplot(pg, aes(x=asm_total_len_bp/1e6, y=asm_scaffold_N50_Mb, color=kingdom))+
  geom_point(alpha=0.6, size=1)+
  scale_color_manual(values=kingdom_cols)+
  labs(x="Assembly size (Mb)", y="Scaffold N50 (Mb)")+
  theme_nature()+theme(
    legend.position = c(0.98, 0.80),
    legend.justification = c("right", "top"),
    legend.background = element_rect(fill = alpha("white", 0.7),
                                     colour = "grey"),
    legend.key.size = unit(3, "mm"),
    legend.text = element_text(size = 6),
    legend.spacing.y = unit(1, "mm"),legend.title = element_blank()
  )
ggsave("FigS1a.pdf",
       p1B,
       width = 60, height = 60, units = "mm", dpi = 300)
## =========================
## 1C: BUSCO completeness
## =========================

stat_busco <- pg %>%
  wilcox_test(BUSCO_completeness ~ kingdom, p.adjust.method = "BH") %>%
  mutate(label = p_label(p.adj))

ymax_busco <- max(pg$BUSCO_completeness, na.rm = TRUE)
stat_busco <- stat_busco %>%
  add_y_position(fun = "max", step.increase = 0.05) %>%  
  mutate(y.position = ymax_busco * c(1.04, 1.08, 1.12))    


p1C <- ggplot(pg, aes(x = kingdom, y = BUSCO_completeness, fill = kingdom, color = kingdom)) +
  geom_jitter(width = 0.15, alpha = 0.6, size = 1.0) +
  geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.75) +
  scale_fill_manual(values = kingdom_cols) +
  scale_color_manual(values = kingdom_cols) +
 # scale_y_continuous(limits = c(60,100)) +
  labs(x = NULL, y = "BUSCO completeness (%)") +
  stat_pvalue_manual(
    stat_busco,
    label = "label",
    y.position = "y.position",
    tip.length = 0.01,
    bracket.size = 0.3,
    size = 2.0,
    inherit.aes = FALSE
  ) +
  theme_nature() 

## =========================
## 1D: TE proportion
## =========================
stat_te <- pg %>%
  wilcox_test(TE_total_pct ~ kingdom, p.adjust.method = "BH") %>%
  mutate(label = p_label(p.adj))

ymax_te <- max(pg$TE_total_pct, na.rm = TRUE)
stat_te <- stat_te %>%
  add_y_position(fun = "max", step.increase = 0.05) %>% 
  mutate(y.position = ymax_te * c(1.10, 1.18, 1.26))

p1D <- ggplot(pg, aes(x = kingdom, y = TE_total_pct, fill = kingdom, color = kingdom)) +
  geom_jitter(width = 0.15, alpha = 0.6, size = 1.0) +
  geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.75) +
  scale_fill_manual(values = kingdom_cols) +
  scale_color_manual(values = kingdom_cols) +
 # scale_y_continuous(expand = expansion(mult = c(0.05, 0.35))) +
  labs(x = NULL, y = "TE proportion (%)") +
  stat_pvalue_manual(
    stat_te,
    label = "label",
    y.position = "y.position",
    tip.length = 0.01,
    bracket.size = 0.3,
    size = 2.0,
    inherit.aes = FALSE
  ) +
  theme_nature()


## =========================
## 1E: Gene density distributions
## =========================

p1E <- ggplot(pg, aes(x=gene_density, y=kingdom, fill=kingdom))+
  geom_density_ridges(alpha=0.8)+
  scale_fill_manual(values=kingdom_cols)+
  labs(x="Gene density (genome)", y=NULL)+
  theme_nature()+theme(
    axis.text.y =element_text(angle = 90, hjust = 0.5, vjust = 0.5))

## =========================
## 1F: hap1–hap2 π distribution
## =========================

p1F <- ggplot(pg, aes(x=pi_genome_wide, color=kingdom))+
  geom_density(linewidth=0.8)+
  scale_color_manual(values=kingdom_cols)+
  labs(x="Genome-wide nucleotide diversity (π)", y="Density")+
  theme_nature()

## =========================
## 1G: hap1–hap2 indel rate
## =========================

p1G <- ggplot(pg, aes(x=indel_rate, color=kingdom))+
  geom_density(linewidth=0.8)+
  scale_color_manual(values=kingdom_cols)+
  labs(x="Genome-wide indel rate", y="Density")+
  theme_nature() + guides(color = guide_legend(title = NULL)) +
  theme(
    legend.position = c(0.55, 0.9),      # 图内右上（可微调）
    legend.justification = c(0, 1),
    legend.direction = "vertical",
    legend.text = element_text(size = 7), # 小一点
    legend.key = element_rect(fill = "transparent", colour = NA),  # 去掉方框
    legend.key.height = unit(3, "mm"),
    legend.key.width  = unit(3, "mm")
  )

## =========================
## combine all panels
## =========================

FigureS1_main <- 
  (p1B | p1C | p1D) /
  ( p1E | p1F | p1G ) +
  plot_annotation(
    tag_levels="a",
    tag_prefix="",
    theme=theme(plot.tag=element_text(face="bold", size=8))
  )

ggsave("FigureS1_main.png", FigureS1_main,
       width=180, height=120, units="mm", dpi=300)
ggsave("FigureS1_main.pdf", FigureS1_main,
       width=180, height=120, units="mm", dpi=300)
