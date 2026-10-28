# Data/inputs

This directory is **intentionally untracked**. Populate it locally from the public accession described in the top-level `README.md`.

## Active input groups

- **SSDO compilation inputs:** curated SSDO / disease-orthologue source files under `ssdo_inputs/` plus the local orthologue tables.
- **Main metric and SA inputs:** `SBGEandSSSdataForMBE.csv`, `FileS7.txt`, `Assoc_Data.txt`, `Ruzicka_et_al_pbio.3000244.s008.csv`, `Innocenti_and_Morrow_ANTAGONISTIC.csv`.
- **Grieshop et al. (2025) FigS3/TableS2 candidates:** copy `DE.summary.tsv` here with its original `FlyBaseID` and `Sig` columns. The shared loader selects the 695 unique FBgn IDs with `Sig == TRUE`; plotting/testing then retains protein-coding genes with mapped Ruzicka estimates.
- **Mishra / SBGE / ASE inputs:** `DifferentialGeneExpression.whole.bodies.tsv`, `Sex.Reversed.AI.test.tsv`, `Heads.FDR.test.for.AI.tsv`, `Gonads.FDR.test.for.AI.tsv`, `Genes.with.Sex.Reversed.AI.tsv`, `Genes.for.analysis.final.tsv`.
- **Fig4 auxiliary inputs:** `dimorphic.subset.list.txt`, `all_paralogs_correct.csv`.
- **Reference annotation:** `Drosophila_melanogaster.BDGP6.54.115.gtf`.

Keep the filenames unchanged so the active scripts can discover them without further editing.
