# Data wrangling

This directory contains the **active build scripts** that create canonical intermediate tables used by the retained plotting and analysis workflow.

## Active scripts

1. `compile_ssdo_reference_data.R`  
   Compiles the canonical SSDO, disease-orthologue, orthologue, and gene-reference tables.

2. `build_mishra_whole_body_sbge_bins.R`  
   Rebuilds the five-bin whole-body SBGE classification used in the current manuscript figures.

3. `build_affx_dros_2_array_reference.R`  
   Builds the Affymetrix Drosophila Genome 2.0 probe reference used to define the Innocenti and Morrow (2010) assay universe for FigS1.

4. `build_sa_resolution_inputs.R`  
   Assembles the study-matched Fig4 feature-membership table, feature-set table, and feature-universe table. FigS7 also uses the membership table to define Penultimate (MB or TS) and Ultimate (HMB or HTS) stages. This is the most involved wrangling step in the repository.

## Brief notes on the more involved transformations

`load_grieshop_candidate_ids.R` is a shared input helper sourced by FigS3 and TableS2, not a standalone build step. It validates `Data/inputs/DE.summary.tsv` and selects the 695 unique significant Grieshop et al. (2025) genes using `Sig == TRUE`.

- **SSDO compilation:** the SSDO, disease-orthologue, and orthologue sets are standardized onto local FBgn-based reference tables so downstream scripts use a single canonical gene-ID space.
- **Affymetrix universe:** the Innocenti and Morrow enrichment tests use the genes represented on the Affymetrix Drosophila Genome 2.0 array rather than the broader protein-coding background.
- **Fig4 feature assembly:** `build_sa_resolution_inputs.R` combines Mishra-derived whole-body SBGE classes, tau-quintile classes (`HTS`, `TS`, `TG`, `HTG`), ASE classes, dimorphic splicing, and duplication into one study-matched feature-membership table. Mechanism columns are `TRUE`/`FALSE` within the corresponding source-study universe and `NA` outside it. The SSS workflow now uses the locally supplied files directly: `Data/inputs/dimorphic.subset.list_wDsRedgenes.txt` for the final DsRed-restored positive set, `Data/inputs/Mishra.sig.SSS_filt25.txt` for the tested-gene denominator, and `Data/inputs/all.exclude.txt` as an optional consistency check that the tested set minus the excluded genes reproduces the final SSS list. `Dup` remains the only protein-coding-genome fallback.

Outputs from these scripts are written to `Data/outputs/` and then consumed by `Plotting_and_analyses/` or `Sup_plotting_and_analyses/`.
