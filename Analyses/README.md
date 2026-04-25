# Analyses

This directory retains broader analytical scripts that are useful context for the project but are not the canonical path for regenerating the current metric plots and active SA figures.

## Scripts

| Script | Status |
| --- | --- |
| `Final_Script.R` | Legacy broader analysis script; not fully updated to the new layout |
| `data_analysis.R` | Legacy extended analysis script; not fully updated to the new layout |
| `SA probability.R` | Legacy sexually antagonistic probability workflow; not yet rebuilt around the new data layout |
| `sa_candidate_enrichment.R` | Current Fisher enrichment workflow for SSDOs, disease orthologues, and orthologues versus legacy SA candidate gene sets; reports the study-matched universes for Innocenti, Ruzicka, and Neal_Truter_Grieshop |
| `sa_resolution.R` | Current SA-resolution Fisher workflow for SSDOs versus assembled feature sets; uses the protein-coding whole-genome universe, reports the four biased SBGE sub-bins (`HFB`, `FB`, `MB`, `HMB`), Mishra-derived `ASE`, `SS-ASE`, `TS-ASE`, and `DR`, and now exports the main split-background and collapsed-background radar variants while keeping the supplemental candidate plots in their simpler two-group form; `UB` is excluded and `DS` remains omitted from the reported outputs |
