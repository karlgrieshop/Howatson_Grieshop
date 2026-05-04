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
   Assembles the Fig4 feature-membership table and feature-set table. This is the most involved wrangling step in the repository.

## Brief notes on the more involved transformations

- **SSDO compilation:** the SSDO, disease-orthologue, and orthologue sets are standardized onto local FBgn-based reference tables so downstream scripts use a single canonical gene-ID space.
- **Affymetrix universe:** the Innocenti and Morrow enrichment tests use the genes represented on the Affymetrix Drosophila Genome 2.0 array rather than the broader protein-coding background.
- **Fig4 feature assembly:** `build_sa_resolution_inputs.R` combines Mishra-derived whole-body SBGE classes, tau-quintile classes (`HTS`, `TS`, `TG`, `HTG`), ASE classes, dimorphic splicing, and duplication into one feature-membership table. The final SSS input is now the local `Data/inputs/dimorphic.subset.list.txt` file, which includes the DsRed-restored genes.

Outputs from these scripts are written to `Data/outputs/` and then consumed by `Plotting/`, `Sup_plotting/`, or `Analyses/`.
