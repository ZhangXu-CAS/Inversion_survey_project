# Inversion manuscript figure code package (cleaned)

This package reorganizes the figure scripts for the inversion manuscript into:
- `main_figures/`
- `extended_data/`
- `supplementary/`
- shared helpers in `helpers/`

Included:
- Main Figure 2, 3, 4 scripts
- Extended Data Figure 1-5 scripts
- Supplementary Figure 1-6 scripts

Expected input files in the package root (or edit path variables at the top of each script):
- `per_genome.csv`
- `per_inversion.csv`
- `summary_stats_all.csv` (Supplementary Fig. 1)
- `inversion_kaks_shift.csv` (Figure 4 + Extended Data Fig. 5)

Outputs are written under `outputs/`.

Notes:
- Figure 1 main is intentionally excluded because it is Illustrator-drawn.
- Panel names and output names were standardized.
- Supplementary Fig. 6 was kept because it appears in the current robustness workflow.

Verification status:
This package was prepared by static code review and reorganization. The environment used here did not include `R`/`Rscript`, so runtime execution against the real CSV inputs was not possible in this session.
