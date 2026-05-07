# Supplementary plotting

This directory contains the active scripts for the retained **supplementary figure** plotting workflow.

## Active scripts

1. `sbge_metric_plots.R`  
   Writes the SBGE-stratified metric plot-data and permutation-test tables kept in `Data/Summary_tables/SBGE_supplement/`.

2. `figureS2_sa_prob_mag.R`  
   Builds `Sup_figures/FigS2_sa_prob_mag.png` from the active SA probability/magnitude inputs.

3. `figureS3_figureS4_metric_panels.R`  
   Builds:
   - `Sup_figures/FigS3_rmf_tau.png`
   - `Sup_figures/FigS4_rlns_tajD_DoS.png`

4. `figureS5_radar_context.R`  
   Builds:
   - `Sup_figures/FigS5_radar_context.png`, a four-panel SA-resolution radar-context figure:
   - panel A: `SSDOs (Figure 4)`
   - panel B: `Innocenti and Morrow (2010)`
   - panel C: `Ruzicka et al. (2019)`
   - panel D: `Neal, Truter, et al. (unpublished)`
   - `Sup_figures/FigS6_SSDO_candidate_overlay.png`, a single-panel overlay comparing SSDOs directly with the Innocenti, Ruzicka, and Neal/Truter SA-candidate sets

   Panels B-D keep any candidate genes that are also SSDOs in the focal candidate set, and exclude only non-candidate SSDOs from those study-specific comparison universes.

The SBGE-stratified component PNGs are no longer kept as tracked standalone outputs.
