# Final Year Project — Sex-Specific Disease Orthologues in *Drosophila melanogaster*

> [!IMPORTANT]
> **Required local data files are intentionally gitignored and must be copied into `Data/inputs/` before rerunning the full workflow. Do not commit these raw files to Git.**
>
> Use the **Grieshop Lab shared/Data** directory here:
> `/Users/rpr23sxu/Library/CloudStorage/OneDrive-UniversityofEastAnglia/Documents/BIO UEA/Grieshop Lab/Grieshop Lab shared/Data`
>
> Required untracked files:
>
> | Local repo path | Needed for | Shared source |
> | --- | --- | --- |
> | `Data/inputs/SBGEandSSSdataForMBE.csv` | core four-category plots | `.../Data/Singh_and_Agrawal_2023/SBGEandSSSdataForMBE.csv` |
> | `Data/inputs/FileS7.txt` | core four-category plots | source file is `.../Data/Fraïsse_et_al_2019/FileS7.csv` |
> | `Data/inputs/all_orthologs.csv` | `compile_ssdo_reference_data.R` | `.../Data/Grieshop_and_Howatson_2026/Orthologues/all_orthologs.csv` |
> | `Data/inputs/all_disease_orthologs.csv` | `compile_ssdo_reference_data.R` | `.../Data/Grieshop_and_Howatson_2026/Orthologues/all_disease_orthologs.csv` |
> | `Data/inputs/Drosophila_melanogaster.BDGP6.54.115.gtf` | gene/transcript reference and SNP-to-gene mapping | `.../Data/Grieshop_and_Howatson_2026/Dmel_GTF/Drosophila_melanogaster.BDGP6.54.115.gtf` |
> | `Data/inputs/Assoc_Data.txt` | Ruzicka SA probability plot | `.../Data/Ruzicka_et_al_2019/Assoc_Data.txt` |
> | `Data/inputs/true_cmh_results.txt` | Neal_Truter_Grieshop SA probability plot | `.../Data/Grieshop_and_Howatson_2026/Neal_Truter_Grieshop/true_cmh_results.txt` |
> | `Data/inputs/gene_hits_snps_90.csv` | Neal_Truter_Grieshop enrichment test | `.../Data/Grieshop_and_Howatson_2026/Neal_Truter_Grieshop/gene_hits_snps_90.csv` |
> | `Data/inputs/Sex.Biased.Genes.Whole.csv` | Mishra SBGE sources | `.../Data/Mishra_et_al_2024/...` |
>
> **Note:** the repo expects a local `Data/inputs/FileS7.txt`; the shared file I found is `FileS7.csv`, so keep a local converted/renamed copy there and leave it untracked.

Analysis code and curated reference data for comparing sex-specific disease orthologues (SSDOs) against broader disease orthologues, all orthologues, and the *D. melanogaster* background genome.

## Repository structure

```text
Data/
  inputs/        raw local inputs copied from shared lab data
  outputs/       canonical ID tables and compiled derived inputs
  Summary_tables/
  unused/
Data_wrangling/  scripts that compile and harmonize reference datasets
Analyses/        broader or legacy analytical scripts retained for context
Plotting/        active plotting workflow for the core four-category figures
Figures/         output figure PDFs
Sup_plotting/    supplemental SBGE-stratified plotting workflow
Sup_figures/     supplemental SBGE-stratified figure PDFs
```

## Current core workflow

1. Run `Data_wrangling/compile_ssdo_reference_data.R`
2. Run `Plotting/plotting.R`
3. Run `Sup_plotting/sbge_metric_plots.R` for the supplemental Mishra-SBGE-stratified versions of the same five core metrics
4. Run `Analyses/sa_candidate_enrichment.R` for the Fisher enrichment tests
5. Run `Plotting/sa_probability_plots.R` for the current SA probability figures
6. Run `Data_wrangling/build_sa_resolution_inputs.R` and then `Analyses/sa_resolution.R` for the SA-resolution feature assembly, the protein-coding whole-genome Fisher tests, and the two-group radar prevalence figure
7. Run `Data_wrangling/build_affx_dros_2_array_reference.R` if `Affx_Dros_2.0.csv` needs to be rebuilt for the study-matched Innocenti enrichment background

The compile step now builds a standalone SSDO table directly from the sex-specific disease inputs and literature orthology files, rather than inheriting SSDOs from a Mishra SBGE merge.
The plotting step is now a single readable script rather than a thin wrapper around helper files.

## Canonical reference sets used by the core plots

The active four-category workflow reads compiled reference files from `Data/outputs/`:

| File | Role |
| --- | --- |
| `Data/outputs/sexspecific_disease_gene_data.csv` | Compiled SSDO table with FlyBase IDs, gene symbols, disease labels, and source provenance |
| `Data/outputs/sexspecific_disease_ids.csv` | Unique SSDO FlyBase IDs used as the focal class |
| `Data/outputs/all_disease_ortholog_ids.csv` | Canonical disease-orthologue FlyBase IDs, harmonized so SSDOs are included |
| `Data/outputs/all_ortholog_ids.csv` | Canonical orthologue FlyBase IDs, harmonized so disease orthologues are included |
| `Data/outputs/dmel_gene_reference.csv` | Gene-level FlyBase reference parsed from the Ensembl/FlyBase GTF |
| `Data/outputs/dmel_transcript_reference.csv` | Transcript-to-gene key parsed from the same GTF |

## Core metric sources

| Dataset | Use |
| --- | --- |
| `Data/inputs/SBGEandSSSdataForMBE.csv` | Tajima's *D*, `rMF`, and `piNfract` (`rlnd`) |
| `Data/inputs/FileS7.txt` | `DoS` and adult tissue specificity |
| `Data/inputs/Sex.Biased.Genes.Whole.csv` | Mishra whole-body SBGE panel assignments for the supplemental SBGE-stratified figures |
| `Data/inputs/Innocenti_and_Morrow_ANTAGONISTIC.csv` | Innocenti and Morrow SA candidate enrichment input |
| `Data/inputs/Ruzicka_et_al_pbio.3000244.s008.csv` | Ruzicka candidate-gene enrichment input |
| `Data/outputs/Affx_Dros_2.0.csv` | Affymetrix Drosophila 2.0 background for the study-matched Innocenti enrichment analysis |

## Notes

- The default “whole-genome” universe used by the active workflows is now the **protein-coding FBgn gene set**.
- `Sup_plotting/sbge_metric_plots.R` uses the same exclusive four-group definitions as the core plots, then restricts each metric to genes retained after merging to Mishra whole-body SBGE categories.
- Raw source files now live in `Data/inputs/`, compiled reference tables live in `Data/outputs/`, and duplicated or superseded files are kept under `Data/unused/`.
- `Analyses/Final_Script.R`, `Analyses/data_analysis.R`, and `Analyses/SA magnitude.R` are retained as older analysis scripts and were **not** fully modernized in this re-layout.
- `Analyses/sa_candidate_enrichment.R` reports both the whole protein-coding genome Fisher enrichment tests and the study-matched-universe Fisher enrichment analyses for Innocenti, Ruzicka, and Neal_Truter_Grieshop.
- `Analyses/sa_resolution.R` reports the protein-coding whole-genome Fisher family and a two-group radar prevalence figure at `Figures/sa_resolution_radar.png`; the reported SBGE terms are the four biased Mishra sub-bins (`HFB`, `FB`, `MB`, `HMB`), `UB` is omitted, and `DS` is excluded from the reported SA-resolution outputs.
- `Plotting/sa_probability_plots.R` produces the default protein-coding-universe SA probability plots in `Figures/` plus full-universe sensitivity versions in `Sup_figures/`.
- `Data_wrangling/build_affx_dros_2_array_reference.R` rebuilds the real Affymetrix Drosophila 2.0 probe-symbol background in `Data/outputs/` for the narrower Innocenti enrichment analysis.
