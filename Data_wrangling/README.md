# Data wrangling

## Active scripts

- `compile_ssdo_reference_data.R`
- `build_affx_dros_2_array_reference.R`

This script compiles the standalone SSDO table, parses the GTF into gene and transcript reference tables, checks the raw subset relationships among SSDOs, disease orthologues, and all orthologues, and writes the canonical ID files consumed by the core plotting workflow.

`build_affx_dros_2_array_reference.R` builds `Affx_Dros_2.0.csv`, a real Affymetrix Drosophila Genome 2.0 probe-symbol reference used only for the study-matched Innocenti enrichment background.
