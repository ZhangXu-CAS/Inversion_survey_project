#!/usr/bin/env Rscript
root <- normalizePath(dirname(sub('^--file=', '', grep('^--file=', commandArgs(FALSE), value=TRUE)[1])), mustWork=TRUE)
setwd(root)
run <- function(x) { message('\n>>> ', x); status <- system2('Rscript', x); if (status != 0) stop('Failed: ', x, call.=FALSE) }
run('supplementary/Phylogenetic_revision_FINAL_Nature.R')
run('main_figures/Fig2_main_cleaned.R')
run('main_figures/Fig3_main_cleaned.R')
run('main_figures/Fig4_main_cleaned.R')
for (x in sort(list.files('extended_data', pattern='\\.R$', full.names=TRUE))) run(x)
for (x in sort(list.files('supplementary', pattern='\\.R$', full.names=TRUE))) if (!grepl('Phylogenetic_revision_FINAL_Nature.R$', x)) run(x)
