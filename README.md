## This repository provides workflows for characterizing chromosomal inversions using haplotype-resolved assemblies across eukaryotic species.

# Pipeline modules
0 data_prepare/: Data download and preprocessing (4 scripts: .py/.sh).

1 assmbly/: HiFi/Hi-C assembly pipelines (4 scripts: .py/.sh).

2 SV_identify/: SV and inversion identification (12 scripts: .py/.sh).

3 Genome_annotation/: Repeat masking and genome annotation (4 scripts: .py/.sh).

4 Post_analysis/: Post-analysis and evolutionary summaries (13 scripts: .py/.sh).

## Typical usage

Prepare data under 0data_prepare/.

Run assembly workflows under 1assmbly/.

Identify inversions under 2SV_identify/.

Annotate genomes under 3Genome_annotation/.

Perform downstream analyses under 4Post_analysis/.

Plot results using 5plot/.
