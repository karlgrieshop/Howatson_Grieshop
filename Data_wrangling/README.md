# Data wrangling

## Active scripts

- `compile_ssdo_reference_data.R`
- `build_affx_dros_2_array_reference.R`
- `build_sa_resolution_inputs.R`

`compile_ssdo_reference_data.R` reads the raw files from `Data/inputs/`, parses the GTF into gene and transcript reference tables, restricts the active canonical universes to protein-coding FBgn genes, and writes the compiled outputs to `Data/outputs/`.

`build_affx_dros_2_array_reference.R` builds `Data/outputs/Affx_Dros_2.0.csv`, a real Affymetrix Drosophila Genome 2.0 probe-symbol reference used only for the study-matched Innocenti enrichment background.

`build_sa_resolution_inputs.R` assembles the current protein-coding SA-resolution feature table and per-feature gene sets used by `Analyses/sa_resolution.R`, including the four biased Mishra SBGE sub-bins (`HFB`, `FB`, `MB`, `HMB`) but not `UB`, and pulling from `Data/inputs/`, the sibling `MaleLimitedEvo` clone, and the shared Harper and Morrow duplication table where available.
