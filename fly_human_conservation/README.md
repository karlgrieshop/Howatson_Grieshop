# Fly–human protein conservation

This workflow measures protein sequence identity and paired-residue coverage
for the selected human orthologue of each eligible fly gene. The released
alignment and mapping tables are in `alignment_local_20261112_local/` and
`mapping_local_20261112_local/`.

## Inputs and references

The workflow requires the canonical tables in `Data/outputs/` and the two
orthology files in `Data/inputs/`. It downloads fixed references:

- Human GRCh38.p14, Ensembl genebuild 2026-04:
  [Ensembl FTP](https://ftp.ebi.ac.uk/pub/ensemblorganisms/GCA/000/001/405/29/ensembl/2026_04/geneset/)
- *D. melanogaster* BDGP6.54, Ensembl release 115:
  [Ensembl FTP](https://ftp.ensembl.org/pub/release-115/)

The selected human partner is the highest-scoring supplied orthology pair for
each fly gene, with HGNC ID used to break score ties. A fixed HGNC-to-Ensembl
mapping table resolves reference-ID ambiguities. For each species, the
representative protein is the valid Ensembl-canonical protein when unique;
otherwise the longest valid protein is selected, with protein ID breaking ties.

## Rebuild

Run commands from the repository root. Python 3.12 is required.

```bash
bash fly_human_conservation/setup_environment.sh
bash fly_human_conservation/recover_mapping.sh --run-id release
bash fly_human_conservation/align_recovered.sh \
  --mapping-dir fly_human_conservation/mapping_local_release --run-id release
```

`recover_mapping.sh` downloads and verifies the fixed references and applies the
included mapping crosswalk. Each run writes to a new directory and refuses to
overwrite existing results. Reference downloads and individual alignments are
not tracked.

The alignment uses Biopython 1.86 `PairwiseAligner` and its BLOSUM62 matrix.
Alignment is global, with gap-opening penalty 10, gap-extension penalty 0.5,
and terminal gaps penalized. Nonstandard residues are masked for scoring but
retained for coverage. Identity is exact standard-residue matches divided by
all alignment columns; coverage is paired-residue columns divided by the full
length of each protein. See Cock et al. (2009) for the Biopython citation
([doi:10.1093/bioinformatics/btp163](https://doi.org/10.1093/bioinformatics/btp163)).

## Figures and tests

```bash
Rscript Sup_plotting_and_analyses/figureS7_fly_human_conservation.R \
  --alignment-dir fly_human_conservation/alignment_local_release \
  --mapping-dir fly_human_conservation/mapping_local_release
Rscript Sup_plotting_and_analyses/figureS8_conservation_adjusted_sa.R \
  --alignment-dir fly_human_conservation/alignment_local_release
```

Fig. S7 compares means using bootstrap intervals and pairwise permutation tests.
Fig. S8 fits SA-effect models with category alone and with category plus
identity and both coverage measures; inference accounts for shared human
partners. Background genes are excluded because conservation is undefined.

Run the regression tests with:

```bash
fly_human_conservation/.venv_local/bin/python -m unittest discover \
  -s fly_human_conservation -p 'test_*.py' -q
Rscript Sup_plotting_and_analyses/test_conservation_sa_regression.R
```
