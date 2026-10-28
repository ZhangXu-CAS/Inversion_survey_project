# Contrasting regimes of chromosomal inversion evolution across plants and animals

**October 2026 revision**

## Overview

This repository contains  scripts and computational workflows for the standardized detection and comparative analysis of chromosomal inversions between phased haplotypes in plants, animals, and fungi. The workflows integrate established third-party tools for sequencing-data preparation, haplotype-resolved assembly, genome annotation, structural-variant discovery, read-based verification, genomic characterization, statistical analysis and figure generation. **This repository is not a new standalone inversion-calling software package.**


## Repository structure

```text
Inversion_survey_project/
├── README.md                      # This overview
├── LICENSE                        # GPL-3.0
├── 0data_prepare/                 # Public sequence retrieval and preparation
├── 1assmbly/                      # Haplotype assembly and scaffolding (original directory spelling)
├── 2SV_identify/                  # Structural-variant and inversion discovery
├── 3Genome_annotation/            # Repeats and gene annotations
├── 4Post_analysis/                # Genomic and evolutionary post-processing
└── 5Rscript_figures/              # Manuscript analysis tables, trees and figure scripts
    ├── README.md                  # Detailed usage and input reference
    ├── run_all.R                  # Figure-analysis entry-point supplied by authors
    ├── main_figures/
    ├── extended_data/
    ├── supplementary/
    ├── helpers/
    ├── per_genome.csv
    ├── per_inversion.csv
    ├── ...                         # Further inputs described in module README
```

## Quick start: summary data and figure analyses

Clone the repository:

```bash
git clone https://github.com/ZhangXu-CAS/Inversion_survey_project.git
cd Inversion_survey_project/5_Rscript_figures
```

Validate the central summary tables and the archived input files:

```bash
Rscript -e 'g <- read.csv("per_genome.csv"); i <- read.csv("per_inversion.csv"); stopifnot(nrow(g) == 1174L, nrow(i) == 130159L); print(table(g$kingdom)); print(table(i$kingdom))'
```

After installing the R package dependencies listed in [`5Rscript_figures/README.md`](5Rscript_figures/README.md), check the figure entry point and run from its own directory:

```bash
Rscript --version
Rscript run_all.R
```

## End-to-end workflow

| Stage | Directory | Primary role | External inputs/requirements |
|---|---|---|---|
| 0 | `0data_prepare/` | Download and organize publicly archived sequencing runs | Accession metadata, ENA/NCBI SRA, SRA Toolkit |
| 1 | `1assmbly/` | Haplotype-aware assembly, phasing and chromosome-scale scaffolding | PacBio HiFi, Hi-C, hifiasm, YaHS and an HPC environment |
| 2 | `2SV_identify/` | Whole-genome alignment, inversion candidate discovery and filtering | Phased genome assemblies, minimap2, SyRI |
| 3 | `3Genome_annotation/` | Repeats, gene prediction and genomic annotations | Genome FASTA, RepeatModeler2/RepeatMasker, BRAKER3, BUSCO |
| 4 | `4Post_analysis/` | Haplotype divergence, inversion context and coding-sequence summaries | Assemblies, alignments, annotations and inversion catalogue |
| 5 | `5Rscript_figures/` | Statistical summaries, sensitivity analyses and figure generation | Verified tabular inputs, species trees, R packages |


## Reproducibility

- Raw HiFi and Hi-C reads: see the paper's sequencing accession list (Supplementary Table 1) and the originating ENA/NCBI SRA records.
- Published genomes and supplementary input resources: see manuscript Data Availability, associated repositories and provenance tables.
- Manuscript numerical inputs: [`5Rscript_figures/`](5Rscript_figures/), including `per_genome.csv`, `per_inversion.csv` and additional analysis-specific tables.
- Figure-generation details: [`5Rscript_figures/README.md`](5Rscript_figures/README.md).

## Code availability, versioning and licence

The code is distributed under the repository's existing **GNU General Public License v3.0** (see `LICENSE`). Third-party software and data retain their own respective licences and citation requirements.

## Citation

If you use this repository, cite the associated manuscript (bibliographic details and DOI to be added on publication) and the versioned software archive when available. Cite original third-party software and data sources as indicated in the manuscript.
