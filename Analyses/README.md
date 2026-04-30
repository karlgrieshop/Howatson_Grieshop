# Analyses

This directory now contains only the two **active manuscript analysis scripts**.

## Active scripts

1. `sa_candidate_enrichment.R`  
   Runs Fisher’s exact enrichment tests for overlap with candidate sexually antagonistic gene sets and writes:
   - `Sup_figures/FigS1_SA_Candidate_Enrichment.png`
   - `Data/Summary_tables/TableS1_SA_Candidate_Enrichment_Tests.csv`

2. `core_analyses.R`  
   Runs the retained SA-resolution enrichment/permutation workflow and writes:
   - `Figures/Fig4_sa_res_radar.png`
   - `Data/Summary_tables/SA_resolution_enrichment_tests.csv`
   - `Data/Summary_tables/SA_resolution_radar_prevalence.csv`
   - `Data/Summary_tables/SA_resolution_radar_permutation_tests.csv`

Legacy exploratory analysis scripts were moved out of the tracked publication surface and are no longer part of the active rerun workflow.
