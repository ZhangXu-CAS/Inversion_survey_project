# Figure analyses and supporting datasets

This folder contains the **manuscript-associated R figure-analysis layer** for *Contrasting regimes of chromosomal inversion evolution across plants and animals*. 

**Reproducibility status:** Input filenames and directory layout below are taken from the authors' 27 October 2026 code folder. 

## Directory layout

```text
5_Rscript_figures/
├── README.md
├── run_all.R                               # Proposed unified entry point
├── main_figures/                          # Main manuscript figure scripts
├── extended_data/                         # Extended Data figure scripts
├── supplementary/                         # Supplementary Figure scripts
├── helpers/                               # Shared R functions and plotting utilities
├── per_genome.csv                         # Genome-level traits and inversion summaries
├── per_inversion.csv                      # Individual inversion records and statistics
├── inversion_kaks_shift.csv               # Inversion-associated coding-sequence analyses
├── population_summary.xlsx               # Population-level analysis summary
├── junction_kingdom_summary_for_figure.tsv  # Scaffold-junction analysis summary
├── allelic_balance_inversions.tsv         # Inversion-level H1/H2 balance summary
├── allelic_balance_samples.tsv           # Genome/sample-level H1/H2 balance summary
├── plant_mating.txt                       # Plant mating-system reference information
├── plant_mcmctree_sanitized.nwk           # Plant dated phylogeny
├── animal_mcmctree_augmented_verified.nwk # Animal dated phylogeny
└── outputs/                               # Generated files (ignored by Git)
```

## Data dictionary and expected role of input files

| File | Data level | Purpose and expected linkage |
|---|---|---|
| `per_genome.csv` | One genome/sample per row | Genome-wide inversion counts, size summaries, genomic covariates; `sample_id` is the primary linkage key |
| `per_inversion.csv` | One inversion per row | Inversion location, length, relative divergence and genomic context; linked by `sample_id`, with unique `inv_id` |
| `inversion_kaks_shift.csv` | Inversion or gene-derived summary | Coding-sequence statistics associated with inversions; verify join keys and model-specific filtering in script |
| `population_summary.xlsx` | Population-derived summary | Independent population-genetic/ascertainment summaries; check sheet names and fields before reuse |
| `junction_kingdom_summary_for_figure.tsv` | Grouped summary | Scaffold-junction technical-sensitivity plots |
| `allelic_balance_inversions.tsv` | Inversion | Competitive haplotype mapping and inversion-specific allelic balance |
| `allelic_balance_samples.tsv` | Genome/sample | Genome-wide allelic-balance controls |
| `plant_mating.txt` | Species/taxon | Mating-system classification and sensitivity analyses |
| `plant_mcmctree_sanitized.nwk` | Species phylogeny | Plant phylogenetic correction; tip names/branch lengths must match analytical samples |
| `animal_mcmctree_augmented_verified.nwk` | Species phylogeny | Animal phylogenetic correction; includes the expanded animal representation |


## Computing environment

- **Operating system:** Unix-like environment supporting Rscript (e.g. Linux or macOS). Record the actual tested OS/kernel in the release.
- **R:** Install a version compatible with the current scripts. Record the tested version using `Rscript --version` and `sessionInfo()`; **a tested minimum version is not yet specified here**.
- **R package dependencies:** inspect every `.R` file for explicit `library()`, `require()` and `package::function()` calls. The project's earlier R code used packages including `ape`, `phytools`, `phylolm`, `nlme`, `dplyr`, `readr`, `tidyr`, `ggplot2`, `patchwork` and `scales`, **but this historical list must not be treated as a verified dependency list for this folder**.
- **Optional tools:** Python 3 for independent input checks; external tools and HPC resources are needed only for upstream assembly and read-level reproduction (see repository-root README).

Recommended after a successful test: save `renv.lock` (or a locked Conda environment), identify Bioconductor dependencies explicitly, and retain R package `sessionInfo()` for each manuscript figure generation run. The correct package installation procedure must be established from those tested dependencies; do not install guessed versions into a manuscript release.

## Quick start

From the **top-level project repository**:

```bash
cd 5_Rscript_figures
Rscript --version
```

Sanity-check input sizes and taxonomic classes **without altering any data**:

```bash
Rscript -e 'g <- read.csv("per_genome.csv"); i <- read.csv("per_inversion.csv"); stopifnot(nrow(g) == 1174L, nrow(i) == 130159L); print(table(g$kingdom)); print(table(i$kingdom))'
```


List the scripts and inspect the top-level orchestration file:

```bash
find main_figures extended_data supplementary helpers -type f -name '*.R' | sort
sed -n '1,200p' run_all.R
```

Individual figures can be rerun with:

```bash
Rscript main_figures/<verified_figure_script>.R
Rscript extended_data/<verified_extended_data_script>.R
Rscript supplementary/<verified_supplementary_script>.R
```

**Replace each filename placeholder with a real script listed above.** Scripts using `source()` may require launching from `5_Rscript_figures/` rather than an arbitrary working directory. Full sample/HPC runtimes and memory usage are currently **not measured** and should be supplied following a real run.

## Figure-to-code mapping (mandatory before the manuscript release)

Create a `FIGURE_MANIFEST.tsv` with **one row per manuscript panel**, containing the columns:

```text
figure_id    panel    script    input_files    output_files    analysis_method    tested_status
```

For example, a manuscript Fig. 2 correlation panel should identify which script fits the phylogenetic regression, which input CSV/tree it uses, how `β` is computed, which output figure it produces, and whether the reproduced values match the submitted paper. **Do not use outdated `Extended_Data_Fig5`/other filenames to infer the current panel mapping.**

### Phylogenetic analysis specification

The supplied revision scripts use `phytools::phylosig(method="lambda", test=TRUE)` to estimate Pagel's λ, `phytools::phylosig(method="K", test=TRUE, nsim=999)` to estimate Blomberg's K, and `phylolm::phylolm(model="lambda", REML=FALSE)` for the revised PGLS models. The original earlier implementation used `nlme::gls()` and `ape::corPagel()`; it should not be silently substituted for the later models.

For the bivariate analyses corresponding to main manuscript Fig. 2c/2f, both the response and genome-wide π predictor are standardized (z-scores) before PGLS, and reported `β` is the **standardized phylogenetically corrected coefficient**, not the ordinary least-squares plotting-line slope. Multivariate PGLS models use log10 genome size, genome-wide TE content and π as primary predictors; gene density is introduced as a sensitivity covariate in the later historical R script. These specifications must be **verified against the final executable code and input trees** before claiming exact figure replication.


## Full versus lightweight reproduction

**Lightweight:** use the two summary CSVs, the trees, R scripts and additional per-figure inputs to reproduce statistical and graphical outputs on a suitably configured workstation/server.

**Full:** regenerate haplotype assemblies, annotations, inversion calls, HiFi read validation, benchmark simulations, provenance checks and derived summary tables from the original external genomic resources. This requires the full upstream pipelines, accessions, relevant databases, HPC specifications and intermediate data not contained in this figure-only folder. The two levels must not be conflated.


## Reporting issues and citation

Please report problems with operating system, R version, R package versions, command executed, input file hashes and error logs. The associated manuscript and its archived code release should be cited once bibliographic details and a software DOI are available. Code inherits the GPL-3.0 licence from the repository unless a subcomponent specifies otherwise; third-party code/data are subject to their own terms.
