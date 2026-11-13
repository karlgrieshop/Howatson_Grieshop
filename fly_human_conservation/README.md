# Fly-human protein conservation

The current workflow runs locally with Bash and Python, without modules or a
scheduler. It downloads fixed references, maps the existing human partners,
aligns one protein pair per fly gene, and supplies FigS8/FigS9. Run commands
from the repository root. Reference preparation and alignment take a few
minutes after downloads; progress is printed to the terminal. Old mapping,
alignment, and job records are retained unchanged as provenance.

## Reconstruct sequences locally

Use Python 3.12 for the pinned Biopython 1.86 / NumPy 1.26.4 environment:

```bash
bash fly_human_conservation/setup_environment.sh
bash fly_human_conservation/recover_mapping.sh --run-id my_run
bash fly_human_conservation/align_recovered.sh \
  --mapping-dir fly_human_conservation/mapping_local_my_run --run-id my_run
```

Run environment setup only once; it refuses to replace an existing environment.
If Python 3.12 is not on `PATH`, provide its absolute executable path with
`bash fly_human_conservation/setup_environment.sh --python /path/to/python3.12`.
The verified local run used the macOS ARM
[standalone CPython 3.12.16 build](https://github.com/astral-sh/python-build-standalone/releases/download/20261110/cpython-3.12.16%2B20261110-aarch64-apple-darwin-install_only.tar.gz),
SHA-256 `e956e60bc9465d577ef8062059aa05cd3a2665b4b8678a2ed24e432368330b8f`.
Other platforms should use their own Python 3.12 build, not this ARM binary.

`recover_mapping.sh` downloads and verifies the fixed human/fly reference
files and applies only the 628 evidence-supported mappings in the frozen,
approved review. It does not choose new human partners. `align_recovered.sh`
uses the explicitly selected mapping. Existing output directories are never
overwritten: change `my_run` or omit `--run-id` for a UTC timestamp.
The installed environment is ignored and can be recreated; it is not copied
from another machine.

The sequence reconstruction needs the canonical `Data/outputs/` tables plus
the original ignored `Data/inputs/all_orthologs.csv` and
`Data/inputs/all_disease_orthologs.csv`. Downloads require public internet
access. FigS8 can instead be rebuilt directly from retained compact tables,
without those raw inputs, references, or Python packages. FigS9 additionally
requires `Data/inputs/Assoc_Data.txt` and the existing R packages listed in its
script.

For initial preparation without reviewed recoveries, a fresh identity review,
or verified reference reuse:

```bash
bash fly_human_conservation/prepare_references.sh --run-id initial
bash fly_human_conservation/review_mapping_recovery.sh --run-id review \
  --reference-source fly_human_conservation/mapping_local_initial
bash fly_human_conservation/recover_mapping.sh --run-id recovered \
  --reuse-verified-references \
  --reference-source fly_human_conservation/mapping_local_initial
bash fly_human_conservation/align_proteins.sh \
  --mapping-dir fly_human_conservation/mapping_local_initial --run-id initial
```

A new review downloads current HGNC identity records and writes
`recovery_review_local_<run-id>/`; it does not apply its proposals or change
the frozen approved review. Sequence/annotation releases stay fixed.
Reference reuse checks the source manifest, checksums, and input identity.
No elevated privileges or scheduler submission are needed.

## Start here: retained inputs

All paths below are repository-relative and tracked:

| Purpose | File |
|---|---|
| Main analysis: one observation per measured fly gene | `fly_human_conservation/alignment_local_20261112_local/tables/gene_conservation.csv` |
| Eligible, excluded, aligned, and per-metric denominators | `fly_human_conservation/alignment_local_20261112_local/tables/category_counts.csv` |
| Every canonical gene, including mapping exclusions | `fly_human_conservation/mapping_local_20261112_local/tables/gene_mapping_audit.csv` |
| All original orthology pairs, source flags, and selection | `fly_human_conservation/mapping_local_20261112_local/tables/all_pairs.csv` |
| Selected pairs, representatives, recovery, xref and score audits | Other CSVs in `mapping_local_20261112_local/tables/`, except the ignored `protein_candidates.csv` |
| Reference releases/checksums, method, completion markers, output hashes | Each run's `provenance/` directory |
| Canonical categories and protein-coding reference | `Data/outputs/{sexspecific_disease_ids,all_disease_ortholog_ids,all_ortholog_ids,dmel_gene_reference}.csv` |

The gene table includes category, HGNC/ENSG/protein IDs, DIOPT score,
representative-selection rules, one-to-one status, original sequence lengths and
hashes, alignment score, masking counts, raw numerators/denominators, and three
percentage metrics. Preserve these fields rather than remapping/reselecting.
Full-graph one-to-one flags and the workflow's raw subset CSV remain as orthology
provenance only; neither FigS8 nor FigS9 performs one-to-one filtered analyses.

## Local statistics and FigS8

Run from the repository root:

```bash
Rscript Sup_plotting_and_analyses/figureS8_fly_human_conservation.R \
  --alignment-dir fly_human_conservation/alignment_local_20261112_local \
  --mapping-dir fly_human_conservation/mapping_local_20261112_local
```

Both directories are explicit inputs. This command **overwrites** the existing
FigS8 figure and summary-table outputs listed below. The script reads the main
gene table, checks category counts against the mapping attrition, and writes:

- `Sup_figures/FigS8_fly_human_conservation.png`
- `Data/Summary_tables/FigS8_category_coverage.csv`
- `Data/Summary_tables/FigS8_group_bootstrap_summaries.csv`
- `Data/Summary_tables/FigS8_pairwise_permutation_tests.csv`
- `Data/Summary_tables/FigS8_input_manifest.csv`

For whole-alignment identity, fly-protein coverage, and human-protein coverage,
the tables report arithmetic means and 95% percentile bootstrap intervals
(5,000 resamples, reproducible seeds). The plots show observed arithmetic
means, intervals, and measured sample sizes. Each metric has three two-sided
permutation tests of mean differences (10,000 permutations), with BH correction
across the three comparisons within that metric.

The updated FigS8 now uses 9,063 measured genes (558 SSDOs, 4,498 disease
orthologues excluding SSDOs, and 4,007 orthologues excluding disease genes).
The original 657 mapping exclusions were reviewed and 628
recovered; the remaining 29 are listed in `RECOVERY_REVIEW.md`.

### Results

In the recovered main set, disease orthologues excluding SSDOs had higher
whole-alignment identity than orthologues excluding disease genes (35.61% vs
34.75%; difference = +0.86 percentage points, two-sided permutation P =
0.0148, BH-adjusted P = 0.0444). Fly coverage was lower in SSDOs than in
disease orthologues excluding SSDOs (83.82% vs 85.72%; difference = -1.90
points, adjusted P = 0.0141), and higher in disease orthologues excluding
SSDOs than in orthologues excluding disease genes (85.72% vs 85.02%;
difference = +0.70 points, adjusted P = 0.0481). No human-coverage contrast
passed BH adjustment.

The main-set pairwise comparisons share 73, 15, and 45 human gene/protein
partners, respectively. This non-independence may make ordinary gene-level
bootstrap intervals and permutation P-values overconfident. FigS9 accounts for
shared human partners directly through clustered inference. The results show
that measured sequence conservation differs across some groups; they do not
establish that conservation explains differences in other biological
outcomes. Non-significance is not evidence of equivalence.

## FigS9: conservation-adjusted SA comparison

From the repository root, using the completed local alignment directory:

```bash
Rscript Sup_plotting_and_analyses/figureS9_conservation_adjusted_sa.R \
  --alignment-dir fly_human_conservation/alignment_<run-id>
```

FigS9 uses only genes with both a finite per-gene `max(abs(Wald_Stat))` from the
original Ruzicka SNP mapping and measured conservation. Background has no
defined fly-human conservation, so a four-category adjusted model is not
identifiable without inventing values. It is excluded, not assigned zero.
The input `Wald_Stat` is a signed Wald z-statistic (approximately `Effect / SD`).
SA effect is its per-gene maximum absolute value, matching Fig1: absolute values
are taken before the gene maximum, without squaring.

The models fit category alone or category plus whole-alignment identity and
both protein coverage measures (common additive linear slopes). FigS9 plots
estimated means before/after adjustment on identical genes, at the relevant
measured universe's pooled conservation means. Two hypotheses are tested before
and after adjustment: a global effect across all three orthologue categories,
then disease-only versus other orthologues after removing SSDOs. BH adjustment
covers these four tests.

Shared human partners and unequal residual variances are handled by a
restricted wild-cluster bootstrap: randomly sign-flip null-model residuals
within human-partner clusters, refit, and compare studentized cluster-robust
Wald statistics (10,000 resamples, plus-one P). This is an approximate robust
regression test, not exact unrestricted residual/label permutation. Mean and
contrast intervals use 5,000 cluster-bootstrap refits. The joint regression
avoids the invalid shortcut of treating pooled response residuals as an
ordinary category outcome without also adjusting category for the covariates.
See `Sup_plotting_and_analyses/README.md` for assumptions, tests, and exports.

The current SA/conservation universe is 8,375/9,063 aligned genes: 530/558
SSDOs, 4,152/4,498 disease-only orthologues, and 3,693/4,007 other orthologues;
the remaining aligned genes lack mapped SA estimates.

On the main matched universe, the disease-only minus other-orthologues
difference is +0.124 before adjustment (95% cluster-bootstrap CI 0.089 to
0.159) and +0.128 after adjustment (0.094 to 0.163); both tests and both
global tests have BH-adjusted P = 0.000100 (the Monte Carlo resolution
floor). The contrast is not attenuated by these measured linear conservation
covariates.

These are observational associations, not causal mediation results. Nonlinear
conservation effects, SNP coverage effects on a maximum statistic, and
gene-family structure beyond the retained human-partner clusters are not
eliminated by this model. The exported `FigS9_*` audits retain the gene-level
SA effect, SNP counts, included/excluded genes, fitted values, residuals,
model coefficients, and exact input/script checksums.

## Measured universes and mapping exclusions

The recovery review completed on 2026-11-11. Mapping job 4955446 applied 628
evidence-supported recoveries; alignment job 4955609 measured all 9,063
mapping-ready pairs. Twenty-nine cases remain excluded. See
[`RECOVERY_REVIEW.md`](RECOVERY_REVIEW.md) for the per-gene evidence, residual
exclusions, and completed jobs.

| Category | Canonical eligible | Excluded | Aligned/measured |
|---|---:|---:|---:|
| SSDOs | 568 | 10 | 558 |
| Disease excluding SSDOs | 4,512 | 14 | 4,498 |
| Orthologues excluding disease | 4,012 | 5 | 4,007 |

These denominators are identical across the three metrics, with no additional
alignment dropout. The R script requires the matching mapping directory
and checks its eligible, excluded, and ready counts against the alignment table.

Of the original 657 pre-alignment exclusions, 628 received reference-aware
mapping recovery. Twenty-nine remain excluded: 18 placeholder `HGNC:0`
identifiers, nine SSDOs absent from the primary pair tables, and two selected
human loci classified as non-protein-coding. Three of the nine have candidate
pairs in ignored local disease-specific orthology CSVs; they were not on the
HPC clone and have not been adjudicated or applied. See `RECOVERY_REVIEW.md`
for the records and follow-up constraints. The 906 LRG xref rows remain a
separate namespace, not additional direct ENSG candidates.

## Fixed scientific definitions

- Categories reproduce tracked canonical IDs: SSDOs; disease IDs minus SSDOs;
  orthologue IDs minus all disease IDs. Background has no retained orthologue
  and is not part of this comparison.
- From the complete supplied pair table, choose the highest DIOPT score per
  fly gene, breaking equal scores by lexicographically sorted stable HGNC ID.
  Do not restrict disease-category genes to disease-annotated human partners,
  use sequence similarity for selection, or rescue a failed winner using a
  lower-scoring pair. The disease input's 10,574 unique pairs are a
  score-consistent subset of the full 35,107 pairs (35,047 in the canonical
  universe). Observed scores span 3-14; no additional threshold was applied,
  and this does not establish upstream high-confidence filtering.
- Exact direct HGNC Symbol xrefs map to ENSG IDs. For the 628 recovered cases,
  the unchanged selected HGNC partner was resolved using authoritative HGNC
  identity records plus corroborating fixed-release annotation, direct xrefs,
  and protein validation. The 18 `HGNC:0` placeholders remain unresolved. LRG
  records are separately audited. Validate the gene/transcript/protein chain
  against GTF and FASTA.
- For each species, choose a uniquely mapped valid Ensembl-canonical protein;
  otherwise choose the longest valid protein, breaking ties by stable protein
  ID. One fly and two human proteins in the retained set used this fallback.
  Strip one terminal stop, retain versioned IDs, and exclude internal stops,
  conflicting mappings, and duplicate FASTA stable IDs.
- Exact global Gotoh alignment: Biopython 1.86, NumPy 1.26.4, BLOSUM62, affine
  gap cost `10 + (gap_length - 1) * 0.5`, including terminal gaps. Fly is the
  target; take the first optimal traceback from the pinned implementation,
  without choosing among ties by conservation. The exact matrix and parameters
  are in alignment provenance.
- Mask positions outside the 20 standard amino acids to X **for scoring only**.
  Original residues remain in saved alignments; no masked position contributes
  an identity match. Masking affects 308 pairs: 33 SSDOs, 158 disease-only,
  117 orthologue-only. It does not exclude these pairs.
- `identity_percent = 100 * standard_identical_columns / alignment_columns`.
  The denominator includes gaps and masked positions: this is conservative
  whole-alignment identity, not BLAST HSP identity.
- Fly/human coverage is `100 * residue_paired_columns / respective length_aa`.
  Masked positions count as residues, not gaps. Coverage is paired-residue
  coverage, not the span covered by a global alignment (which would be 100%).
- One-to-one means degree one at both ends of the **full supplied bipartite
  orthology graph**, before canonical filtering, mapping, or alignment. Keep
  the existing flag; do not redefine it after exclusions.

## References, jobs, and provenance

Human: GRCh38.p14 / GCA_000001405.29 / Ensembl genebuild 2026-04, using
`pep.fa.bgz`, `genes.gtf.bgz`, and `xref.tsv.gz` from
<https://ftp.ebi.ac.uk/pub/ensemblorganisms/GCA/000/001/405/29/ensembl/2026_04/geneset/>.
Fly: BDGP6.54 / Ensembl release 115, using
`Drosophila_melanogaster.BDGP6.54.pep.all.fa.gz` and
`Drosophila_melanogaster.BDGP6.54.115.gtf.gz`. URLs, retrieval dates, provider
MD5/BSD checks, and SHA-256 fingerprints are retained in reference manifests
and `prepare_references.py`.

| Historical job | Outcome |
|---|---|
| 4935248 preparation | Failed on an LRG namespace assumption, after all five reference files passed checksums; original scripts/logs retained |
| 4935338 preparation recovery | Completed, exit 0, 77 seconds; reused/reverified references; fresh mapping outputs |
| 4935733 global alignment | Completed, exit 0, 3m57s, about 1.36 GiB peak RSS; all 8,435 pairs aligned |
| 4955170 recovery review | Completed, exit 0; identified 628 supported proposals and 29 remaining exclusions |
| 4955385 recovery mapping | Failed on no-pair flag validation; partial outputs and logs retained |
| 4955446 recovery mapping | Completed, exit 0; applied 628 approved mappings |
| 4955609 recovered alignment | Completed, exit 0, 3m48s, about 1.36 GiB peak RSS; all 9,063 pairs aligned |

On 2026-11-12 the entire recovered sequence workflow ran locally from freshly
downloaded fixed references, using CPython 3.12.16, Biopython 1.86, and NumPy
1.26.4 on macOS ARM. `mapping_local_20261112_local/` took 123 seconds and
`alignment_local_20261112_local/` took 101 seconds. All 9,063 alignment rows
and biological mapping tables matched the recovered historical outputs
exactly. FigS8 and FigS9 were then regenerated from these local outputs.

A fresh full local review, `recovery_review_local_20261112_local_verify/`,
independently reproduced all 657 review decisions and proposed IDs: 628
supported recoveries, 18 identifier-history cases, nine no-pair cases, and two
classified non-protein-coding exclusions. Its newly retrieved HGNC snapshot
had the same SHA-256 as the approved snapshot. It was not silently applied as
a new decision rule. See `RECOVERY_REVIEW.md` for the additional source-pair
leads that require a separate orthology-universe decision.

Historical CSVs, JSONs, logs, completion markers, and manifests were moved
**without rewriting bytes**. They intentionally record the original
`side/fly_human_conservation/` paths. Resolve that prefix to
`fly_human_conservation/` for the relocated files; `Data/` paths did not change.
Relative alignment paths inside each gene table are relative to that job's
alignment directory and point to ignored per-pair files.
The directory's `.gitattributes` preserves original CSV/TSV line endings,
empty fields, and archived provenance bytes during checkout.

The original executed source versions are in `provenance/executed_scripts/`;
their hashes match the successful jobs' provenance. The failed version is
in `failed_job_4935248_scripts/`. These are **archives, not runnable entry
points**. `provenance/HPC_EXECUTION_NOTES.md` preserves the chronological
historical notes, with their original paths. Current entry scripts run locally;
historical source snapshots retain their original scheduler-specific code.

## Tracked versus local bulk files

Git includes scripts/tests/requirements, compact tables, reference and output
manifests, completion records, original source snapshots, and successful/failed
logs. It excludes:

- Reference-download directories under `fly_human_conservation/`: compressed
  references and provider checksum files, including run-specific downloads.
- `fly_human_conservation/alignment_*/alignments/`: individual sequence/traceback
  JSONs (listed by checksum in the tracked output manifest).
- `fly_human_conservation/mapping_*/tables/protein_candidates.csv`: large
  protein-candidate audit, not needed for statistics on the approved gene set.
- `fly_human_conservation/.venv*/`: installed local Python environment and
  optional standalone interpreter. The old ignored environment under `side/`
  is no longer used.
- The original `Data/inputs/all_orthologs.csv` and
  `Data/inputs/all_disease_orthologs.csv`, still ignored. The tracked full-pair
  table already preserves their pair content and source membership.

Compact tables are sufficient for FigS8; FigS9 also reads the original SA
SNP data. Sequence reconstruction needs the raw orthology inputs and reference
downloads. Never copy installed environments between machines; recreate them
from `alignment_requirements.txt`.

## Verification

```bash
mkdir -p fly_human_conservation/.venv_local/scratch
TMPDIR="$PWD/fly_human_conservation/.venv_local/scratch" PYTHONDONTWRITEBYTECODE=1 \
  fly_human_conservation/.venv_local/bin/python -m unittest discover \
  -s fly_human_conservation -p 'test_*.py' -q
Rscript Sup_plotting_and_analyses/test_conservation_sa_regression.R
```

The local Python suite has 22 tests. Preparation and alignment completion
records identify their output-manifest hashes; manifests retain hashes for
compact tables and ignored bulk outputs. `check_handoff.py` remains a
standard-library check of the **original** archived handoff, not a substitute
for validating a fresh run's manifests. Terminal progress can be redirected
to a user-chosen log file; per-run execution details are always retained under
`provenance/`.
