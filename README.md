# Final Year Project — Sex-Specific Disease Orthologues in *Drosophila melanogaster*

> [!IMPORTANT]
> **Required local data files are intentionally gitignored and must be copied into `Data/` before rerunning the full workflow. Do not commit these raw files to Git.**
>
> Use the **Grieshop Lab shared/Data** directory here:
> `/Users/rpr23sxu/Library/CloudStorage/OneDrive-UniversityofEastAnglia/Documents/BIO UEA/Grieshop Lab/Grieshop Lab shared/Data`
>
> Required untracked files:
>
> | Local repo path | Needed for | Shared source |
> | --- | --- | --- |
> | `Data/SBGEandSSSdataForMBE.csv` | core four-category plots | `.../Data/Singh_and_Agrawal_2023/SBGEandSSSdataForMBE.csv` |
> | `Data/FileS7.txt` | core four-category plots | source file is `.../Data/Fraïsse_et_al_2019/FileS7.csv` |
> | `Data/all_orthologs.csv` | `compile_ssdo_reference_data.R` | `.../Data/Grieshop_and_Howatson_2026/Orthologues/all_orthologs.csv` |
> | `Data/all_disease_orthologs.csv` | `compile_ssdo_reference_data.R` | `.../Data/Grieshop_and_Howatson_2026/Orthologues/all_disease_orthologs.csv` |
> | `Data/Drosophila_melanogaster.BDGP6.54.115.gtf` | gene/transcript reference, background gene universe, SNP-to-gene mapping | `.../Data/Grieshop_and_Howatson_2026/Dmel_GTF/Drosophila_melanogaster.BDGP6.54.115.gtf` |
> | `Data/Assoc_Data.txt` | Ruzicka SA probability plot | `.../Data/Ruzicka_et_al_2019/Assoc_Data.txt` |
> | `Data/true_cmh_results.txt` | Neal_Truter_Grieshop SA probability plot | `.../Data/Grieshop_and_Howatson_2026/Neal_Truter_Grieshop/true_cmh_results.txt` |
> | `Data/gene_hits_snps_90.csv` | Neal_Truter_Grieshop enrichment test | `.../Data/Grieshop_and_Howatson_2026/Neal_Truter_Grieshop/gene_hits_snps_90.csv` |
> | `Data/Sex.Biased.Genes.Whole.csv` | older legacy analyses in `Analyses/` | `.../Data/Mishra_et_al_2024/...` |
>
> **Note:** the repo expects a local `Data/FileS7.txt`; the shared file I found is `FileS7.csv`, so keep a local converted/renamed copy in `Data/` and leave it untracked.

Analysis code and curated reference data for comparing sex-specific disease orthologues (SSDOs) against broader disease orthologues, all orthologues, and the *D. melanogaster* background genome.

## Repository structure

```text
Data/            curated inputs, canonical ID files, and summary tables
Data_wrangling/  scripts that compile and harmonize reference datasets
Analyses/        broader or legacy analytical scripts retained for context
Plotting/        active plotting workflow for the core four-category figures
Figures/         output figure PDFs
```

## Current core workflow

1. Run `Data_wrangling/compile_ssdo_reference_data.R`
2. Run `Plotting/plotting.R`
3. Run `Analyses/sa_candidate_enrichment.R` for the Fisher enrichment tests
4. Run `Plotting/sa_probability_plots.R` for the current SA probability figures

The compile step now builds a standalone SSDO table directly from the sex-specific disease inputs and literature orthology files, rather than inheriting SSDOs from a Mishra SBGE merge.
The plotting step is now a single readable script rather than a thin wrapper around helper files.

## Canonical reference sets used by the core plots

The active four-category workflow reads these files from `Data/`:

| File | Role |
| --- | --- |
| `sexspecific_disease_gene_data.csv` | Compiled SSDO table with FlyBase IDs, gene symbols, disease labels, and source provenance |
| `sexspecific_disease_ids.csv` | Unique SSDO FlyBase IDs used as the focal class |
| `all_disease_ortholog_ids.csv` | Canonical disease-orthologue FlyBase IDs, harmonized so SSDOs are included |
| `all_ortholog_ids.csv` | Canonical orthologue FlyBase IDs, harmonized so disease orthologues are included |
| `dmel_gene_reference.csv` | Gene-level FlyBase reference parsed from the Ensembl/FlyBase GTF |
| `dmel_transcript_reference.csv` | Transcript-to-gene key parsed from the same GTF |

## Core metric sources

| Dataset | Use |
| --- | --- |
| `Data/SBGEandSSSdataForMBE.csv` | Tajima's *D*, `rMF`, and `piNfract` (`rlnd`) |
| `Data/FileS7.txt` | `DoS` and adult tissue specificity |
| `Data/Drosophila_melanogaster.BDGP6.54.115.gtf` | Reference gene and transcript identifiers for background and ID harmonization |

## Notes

- The core plots now maximise inclusion by assigning SSDOs from the compiled SSDO table, not from a pre-merged SBGE subset.
- Raw source files that are duplicated, superseded, or not part of the active four-category workflow are kept under `Data/unused/`.
- `Analyses/Final_Script.R`, `Analyses/data_analysis.R`, and `Analyses/SA magnitude.R` are retained as older analysis scripts and were **not** fully modernized in this re-layout.
- `Analyses/sa_candidate_enrichment.R` restores the Fisher enrichment tests against legacy SA candidate gene sets using the current canonical SSDO IDs.
- `Plotting/sa_probability_plots.R` restores the SA probability plots from the original study-specific inputs now present in `Data/`, using gene-level maxima after mapping SNPs onto the local gene reference.
