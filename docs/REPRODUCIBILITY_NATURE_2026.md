# Nature revision: code and reproducibility status

**Status: documentation and lightweight audit setup; NOT a verified full-paper reproduction.**
This branch does **not** replace the canonical `main` datasets, regenerate manuscripts, or certify existing results.

## Manuscript
*Contrasting regimes of chromosomal inversion evolution across plants and animals* (October 2026 revision).

## Authoritative inputs
The final sample-level and inversion-level files are named `per_genome.csv` and `per_inversion.csv`. The **author-confirmed** files (uploaded 2026-10-23) have:
- per_genome: 1,174 genomes = 982 Metazoa + 181 Viridiplantae + 11 Fungi.
- per_inversion: 130,159 inversion rows = 119,022 Metazoa + 10,976 Viridiplantae + 161 Fungi.

**Important:** The `5figure_Rplot/per_genome.csv` and `5figure_Rplot/per_inversion.csv` currently on `main` are older files of different sizes. They must not be presented as the manuscript-final datasets until authors verify replacement hashes. Upload and verify authoritative files before attempting final reproduction.

## Lightweight data audit (no HPC)
Requires **Python 3.9+**; uses the standard library. Copy the two verified final CSVs into a directory, then:

```bash
python scripts/audit_final_tables.py --genome /path/to/per_genome.csv --inversion /path/to/per_inversion.csv --json audit.json
```

Exit 0 indicates that the declared sample and inversion counts, identifiers and sample-level inversion totals are consistent. Exit 1 indicates failing checks; the report includes SHA-256 hashes and warnings. **Passing does not validate the scientific inference, primary BAMs, or figures.**

A small deterministic fixture can be run in CI with `python -m unittest discover -s tests -v`. Typical runtime on the full 130k-inversion CSV must be measured on the deployment machine; no unmeasured installation/runtime claim is made.

## Analysis pipelines and known upstream dependencies

| Component | Current source/implementation | Missing for final reproduction |
|---|---|---|
| Acquisition and assembly | `0data_prepare`, `1assmbly` | ENA/SRA accessions, full dependencies, HPC recipe |
| SV calling | `2SV_identify` | Final parameter manifest, read validation evidence |
| Genome annotation and evolution | `3Genome_annotation`, `4Post_analysis` | Verified input references and output hashes |
| Main/Extended/Supplementary figures | `5figure_Rplot` | Final October figure numbering and verified scripts |
| Phylogenetic signal | Historical `Phylogenetic_revision_FINAL_Nature.R` | Audited R script + final dated trees + sessionInfo |
| PGLS | Historical `Phylogenetic_revision_FINAL_Nature.R` | Refit to 1,174 genomes and release results |
| Simulation, provenance, heterozygote ascertainment | Methods/Supplementary Notes | Full scripts, seed, parameter logs and source inputs |

### Phylogenetic models, from preserved R scripts

The more recent historical R script calls:
- `phytools::phylosig(tree, trait, method="lambda", test=TRUE)` for Pagel's λ;
- `phytools::phylosig(tree, trait, method="K", test=TRUE, nsim=999)` for Blomberg's K;
- `phylolm::phylolm(..., model="lambda", REML=FALSE)` for bivariate and multivariate PGLS.
- `ape` for matching/pruning phylogenies and `stats::p.adjust(..., method="BH")` for adjustment.

**Fig. 2c and 2f:** Response and genome-wide π predictor are separately standardized to z-scores, and the reported β is the standardized PGLS regression coefficient. Ordinary least-squares visualization lines are not identical to PGLS effects.

**Multivariate PGLS:** the primary covariates are log10 genome size, TE fraction and genome-wide π, with gene density added as a **sensitivity** covariate. Plant and animal trees are analysed separately.

**Do not reuse the historical script unmodified:** its hardcoded initial assertions require 1,180 total / 187 plant genomes, which differ from authoritative 1,174 / 181. Final results, including confidence intervals and tree-tip mapping, need to be recalculated and compared to manuscript values. The historical `Phylogenetic_signal_revision.R` instead used `nlme::gls` with `ape::corPagel` and is not interchangeable with the final `phylolm` workflow.

## Verified vs unverified
- **Static inspection only:** code location and legacy functions; the repository currently contains both older figure numbering and previous datasets.
- **Not yet verified:** executing R PGLS or full assembly, reproducing Nature figures, HiFi/BAM-level validation, dated-tree branch scales, HPC performance, missingness sensitivity, final reproducibility claims.
- **Do not hide:** 71 records with inversion length <200 bp; mixed log definitions in former mean-length results; legacy scaffold/source counts; numerous exactly-one flank ratios.

## Nature code/software checklist
Custom Bash/Python/R workflows should be declared as custom analysis code, not new inversion-calling software. License in the repo is GPL-3.0. The public repository address is https://github.com/ZhangXu-CAS/Inversion_survey_project . A **version-specific DOI cannot be reported until a tested release is archived**. Referees must be supplied with code and usage directions; any unverified steps must be identified as such.

## Before archival release
1. Replace older public CSVs with author-approved byte-identical final CSVs, recording SHA-256.
2. Add the actual final phylogenetic script and dated trees, and run them in R with logged package versions.
3. Map every manuscript panel to versioned input, script, output and test.
4. Add demo data with expected outputs; run and document both a smoke test and an HPC command.
5. Record final code commit SHA; create a release and DOI-backed snapshot; update Code Availability only after verification.
