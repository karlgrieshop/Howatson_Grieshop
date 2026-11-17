# SSDO evolutionary genetics analyses

This repository contains analysis scripts and retained results for the study of
sex-specific disease orthologues (SSDOs) in *Drosophila melanogaster*. Analyses
compare SSDOs, other disease orthologues, other orthologues, and background
genes.

## Repository contents

| Path | Contents |
|---|---|
| `Data/inputs/` | Local input files; excluded from Git |
| `Data/outputs/` | Compiled gene sets and reference tables |
| `Data/Summary_tables/` | Statistical summaries and figure data |
| `Data_wrangling/` | Scripts that build analysis inputs |
| `Plotting_and_analyses/` | Main-text analyses and figures |
| `Sup_plotting_and_analyses/` | Supplementary analyses and figures |
| `Figures/`, `Sup_figures/` | Manuscript figures |
| `Tables/` | Table-generation scripts and editable tables |
| `fly_human_conservation/` | Fly–human sequence workflow and results |

## Reproduction

Run R scripts from the repository root. Required input filenames and source
studies are listed in [`Data/inputs/README.md`](Data/inputs/README.md).
Directory-specific READMEs describe each analysis and its outputs.

The fly–human sequence workflow uses Python 3.12 with the pinned packages in
`fly_human_conservation/alignment_requirements.txt`. Its commands, reference
releases, methods, and verification steps are in
[`fly_human_conservation/README.md`](fly_human_conservation/README.md).

Tracked summary tables and figures are the released analysis outputs. Running
the scripts may overwrite the corresponding outputs.
