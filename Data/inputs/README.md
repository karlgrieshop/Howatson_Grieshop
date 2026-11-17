# Input data

Input files are excluded from Git and available at
[10.5281/zenodo.20162756](https://doi.org/10.5281/zenodo.20162756).
Extract `input_data.zip` into the repository root to populate this directory.
Input filenames and source studies are listed below.

## Gene sets and references

- `ssdo_inputs/sex_specific_diseases/*.txt`
- `ssdo_inputs/orthology_sources/diopt_results_PCOS_genes.csv`
- `ssdo_inputs/orthology_sources/diopt_results.csv`
- `ssdo_inputs/orthology_sources/coronary_artery_disease_orthologs.csv`
- `ssdo_inputs/orthology_sources/osteoporosis_orthologs.csv`
- `ssdo_inputs/orthology_sources/alzheimer_disease_orthologs.csv`
- `ssdo_inputs/orthology_sources/parkinsons_orthologs.csv`
- `all_disease_orthologs.csv`, `all_orthologs.csv`
- `Drosophila_melanogaster.BDGP6.54.115.gtf`

## Expression, population-genetic, and association data

- `SBGEandSSSdataForMBE.csv` (Singh and Agrawal, 2023) (doi:10.1093/molbev/msad091)
- `FileS7.txt` (Fraïsse et al., 2019) (doi:10.1093/molbev/msy246)
- `Assoc_Data.txt`, `Ruzicka_et_al_pbio.3000244.s008.csv` (Ruzicka et al., 2019) (doi:10.1371/journal.pbio.3000244)
- `Innocenti_and_Morrow_ANTAGONISTIC.csv` (Innocenti and Morrow, 2010) (doi:10.1371/journal.pbio.1000335)
- `DE.summary.tsv` (Grieshop et al., 2025) (doi:10.1093/molbev/msaf099)
- `dimorphic.subset.list_wDsRedgenes.txt` `all.exclude.txt` (modified from Grieshop et al., 2025) (doi:10.1093/molbev/msaf099)
- `DifferentialGeneExpression.whole.bodies.tsv`, `Sex.Reversed.AI.test.tsv`, `Heads.FDR.test.for.AI.tsv`,
  `Gonads.FDR.test.for.AI.tsv`, `Genes.with.Sex.Reversed.AI.tsv`, `Genes.for.analysis.final.tsv` (Mishra et al., 2024a) (doi:10.1093/gbe/evae234)
- `Mishra.sig.SSS_filt25.txt` (Mishra et al., 2024b) (doi:10.1093/evolut/qpae004)
- `all_paralogs_correct.csv` (Harper and Morrow, 2025) (doi:10.1098/rspb.2024.2629)

Keep the original filenames and file contents. The scripts use these names to
locate inputs.
