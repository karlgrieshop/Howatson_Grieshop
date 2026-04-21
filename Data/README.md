# Data

This directory now separates active reference data from archived or superseded files.

## Active files for the core four-category workflow

| File | Purpose |
| --- | --- |
| `sexspecific_disease_gene_data.csv` | Compiled SSDO table built from sex-specific disease text files and literature orthology sources |
| `sexspecific_disease_ids.csv` | Unique SSDO FlyBase IDs used in the core plots |
| `all_disease_ortholog_ids.csv` | Canonical disease-orthologue ID set |
| `all_ortholog_ids.csv` | Canonical orthologue ID set |
| `dmel_gene_reference.csv` | Gene-level reference parsed from the GTF |
| `dmel_transcript_reference.csv` | Transcript-level reference parsed from the GTF |
| `FileS7.txt` | Fraisse et al. metric source |
| `SBGEandSSSdataForMBE.csv` | Singh & Agrawal metric source |
| `Sex.Biased.Genes.Whole.csv` | Mishra whole-body SBGE table for later stratified analyses |
| `Summary_tables/` | Current plot-ready tables, raw plot data, and permutation-test outputs |

## SSDO input sources

- `ssdo_inputs/sex_specific_diseases/` contains the FlyBase-derived sex-specific disease text files.
- `ssdo_inputs/orthology_sources/` contains the literature and DIOPT-derived orthology source tables used to add SSDOs beyond the FlyBase disease-search lists.

## Archived files

`unused/` holds duplicated text files, superseded intermediate outputs, and legacy data inputs that are not part of the active four-category workflow.
