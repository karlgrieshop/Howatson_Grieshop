# Supplementary plotting

This directory contains the active scripts for the retained **supplementary figure** plotting workflow.

## Active scripts

1. `sbge_metric_plots.R`  
   Writes the SBGE-stratified metric plot-data and permutation-test tables kept in `Data/Summary_tables/SBGE_supplement/`.

2. `figureS2_context.R`  
   Builds `Sup_figures/FigS2_context.png`, a Fig1-style four-panel context figure comparing:
   - SSDOs
   - Innocenti & Morrow candidates
   - Ruzicka et al. candidates
   - Neal, Truter et al. candidates
   
   Each panel also includes horizontal blue/green/grey reference lines marking the original Fig1 bootstrap means for disease orthologues, orthologues, and background.

3. `figureS3_sa_prob_mag.R`  
   Builds `Sup_figures/FigS3_sa_prob_mag.png` from the active SA probability/magnitude inputs.

4. `figureS4_figureS5_metric_panels.R`  
   Builds:
   - `Sup_figures/FigS4_rmf_tau.png`
   - `Sup_figures/FigS5_rlns_tajD_DoS.png`

5. `figureS6_radar_context.R`  
    Builds:
   - `Sup_figures/FigS6_radar_context.png`, a four-panel study-matched SA-resolution radar-context figure:
   - panel A: `SSDOs (Figure 4)`
   - panel B: `Innocenti and Morrow (2010)`
   - panel C: `Ruzicka et al. (2019)`
   - panel D: `Neal, Truter, et al. (unpublished)`
   - `Sup_figures/FigS7_SSDO_candidate_overlay.png`, a study-matched single-panel overlay comparing SSDOs directly with the Innocenti, Ruzicka, and Neal/Truter SA-candidate sets

    Panels B-D keep any candidate genes that are also SSDOs in the focal candidate set, exclude only non-candidate SSDOs from those study-specific comparison universes, and use the same study-matched per-feature denominators as Fig4. In the current tracked workflow, `SSS` uses the coauthor-supplied `Mishra.sig.SSS_filt25.forKarl.txt` tested-gene frame together with the DsRed-restored final list in `dimorphic.subset.list_wDsRedgenes.txt`, and `Dup` remains the only protein-coding-genome fallback.

The SBGE-stratified component PNGs are no longer kept as tracked standalone outputs.
