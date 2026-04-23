# Analyses

This directory retains broader analytical scripts that are useful context for the project but are not the canonical path for regenerating the core four-category figures.

## Scripts

| Script | Status |
| --- | --- |
| `Final_Script.R` | Legacy broader analysis script; not fully updated to the new layout |
| `data_analysis.R` | Legacy extended analysis script; not fully updated to the new layout |
| `SA magnitude.R` | Legacy sexually antagonistic SNP workflow; not yet rebuilt around the new data layout |
| `sa_candidate_enrichment.R` | Current Fisher enrichment workflow for SSDOs versus legacy SA candidate gene sets; uses the whole protein-coding genome as the default broad universe and keeps the parallel study-matched universes |
| `sa_resolution.R` | Current SA-resolution Fisher workflow for SSDOs versus assembled feature sets; uses the protein-coding whole-genome universe, reports the four biased SBGE sub-bins (`HFB`, `FB`, `MB`, `HMB`) rather than aggregate `SBGE`, excludes `UB`, and omits `DS` from the reported outputs |
