# Fly-human conservation: local analysis handoff

**HPC preparation and global alignment are complete.** This remains exploratory,
separate from canonical gene sets and manuscript figures/tables. On 2026-11-10,
the user requested moving the work from ignored `side/fly_human_conservation/`
to this tracked directory. No sequences, metrics, categories, or exclusion rules
were changed during that move.

The current compact inputs are tracked; local statistics and plotting need no
raw sequences, reference downloads, alignments, or Python environment. The
requested mapping recovery requires a fresh HPC output set to be committed or
otherwise made available locally before rerunning FigS8.

From the local repository root, verify the handoff with standard Python only:

```bash
python3 fly_human_conservation/check_handoff.py
```

## Start here: retained inputs

All paths below are repository-relative and tracked:

| Purpose | File |
|---|---|
| Main analysis: one observation per measured fly gene | `fly_human_conservation/alignment_4935733/tables/gene_conservation.csv` |
| Full-graph one-to-one sensitivity subset | `fly_human_conservation/alignment_4935733/tables/one_to_one_conservation.csv` |
| Eligible, excluded, aligned, and per-metric denominators | `fly_human_conservation/alignment_4935733/tables/category_counts.csv` |
| Every canonical gene, including mapping exclusions | `fly_human_conservation/mapping_4935338/tables/gene_mapping_audit.csv` |
| All original orthology pairs, source flags, and selection | `fly_human_conservation/mapping_4935338/tables/all_pairs.csv` |
| Selected pairs, representatives, xref and score audits | Other CSVs in `fly_human_conservation/mapping_4935338/tables/`, except the ignored `protein_candidates.csv` |
| Reference releases/checksums, method, completion markers, output hashes | `fly_human_conservation/{provenance,mapping_4935338/provenance,alignment_4935733/provenance}/` |
| Canonical categories and protein-coding reference | `Data/outputs/{sexspecific_disease_ids,all_disease_ortholog_ids,all_ortholog_ids,dmel_gene_reference}.csv` |

The gene table includes category, HGNC/ENSG/protein IDs, DIOPT score,
representative-selection rules, one-to-one status, original sequence lengths and
hashes, alignment score, masking counts, raw numerators/denominators, and three
percentage metrics. Preserve these fields rather than remapping/reselecting.

## Local statistics and FigS8

Run from the repository root:

```bash
Rscript Sup_plotting_and_analyses/figureS8_fly_human_conservation.R \
  --alignment-dir fly_human_conservation/alignment_4935733
```

The alignment directory is an explicit input so the analysis can be rerun
against a fresh, uniquely named HPC result without overwriting prior outputs.
The script reads its main and full-graph one-to-one gene tables and writes:

- `Sup_figures/FigS8_fly_human_conservation.png`
- `Data/Summary_tables/FigS8_category_coverage.csv`
- `Data/Summary_tables/FigS8_group_bootstrap_summaries.csv`
- `Data/Summary_tables/FigS8_pairwise_permutation_tests.csv`
- `side/fly_human_conservation/one_to_one_sensitivity.png` (exploratory)

For whole-alignment identity, fly-protein coverage, and human-protein coverage,
the tables report arithmetic means and 95% percentile bootstrap intervals
(5,000 resamples, reproducible seeds). The plots show observed arithmetic
means, intervals, and measured sample sizes. Each metric has three two-sided
permutation tests of mean differences (10,000 permutations), with BH correction
across the three comparisons within that metric. The one-to-one sensitivity is
a separate analysis set with its own three-test BH family per metric.

**The currently retained alignment is the initial filtered analysis, not the
final FigS8 dataset.** Mapping recovery is now requested for genes excluded
before alignment. Until the refreshed mapping/alignment outputs are available,
the figure and results below describe only the initial 8,435 measured genes
(527 SSDOs, 4,195 disease orthologues
excluding SSDOs, and 3,713 orthologues excluding disease genes). All three
metrics have these same denominators. The one-to-one subset has 2,299 measured
genes (101, 750, and 1,448, respectively). The 657 mapping exclusions among
9,092 eligible genes were initially omitted, not imputed as zero. The requested
mapping review may recover some; regenerate these counts and all statistics
from the refreshed alignment before treating FigS8 as final.

### Results

In the main analysis, disease orthologues excluding SSDOs had mean identity
35.76%, compared with 34.90% in orthologues excluding disease genes. This
difference was not significant after within-metric BH adjustment
(two-sided permutation P = 0.0217; adjusted P = 0.0651). Mean fly coverage
was 85.66% versus 84.97% in those two groups (P = 0.0446; adjusted P = 0.0549).
The only main-analysis contrast passing BH adjustment was lower fly coverage
for SSDOs (83.58%) than disease orthologues excluding SSDOs (85.66%;
difference = -2.09 percentage points, adjusted P = 0.00750). No human-coverage
contrast passed BH adjustment.

In the one-to-one sensitivity, disease orthologues excluding SSDOs had higher
identity than orthologues excluding disease genes (37.60% versus 35.91%;
difference = 1.69 points, adjusted P = 0.0372) and higher fly coverage
(90.25% versus 87.21%; difference = 3.04 points, adjusted P = 0.000300).
No other one-to-one contrast passed its within-metric BH adjustment; no
human-coverage contrast was significant after adjustment.

The main-set pairwise comparisons share 69, 14, and 41 human gene/protein
partners, respectively. Therefore, fly genes are not fully independent in the
main analysis, and ordinary gene-level bootstrap intervals and gene-label
permutation P-values may be overconfident. The one-to-one sensitivity has no
shared human partners by construction and is an important check, but its SSDO
sample is small (n = 101). These analyses compare average conservation; they
do not establish that conservation explains differences in other outcomes.
Non-significance is not evidence of equivalence.

These are exploratory results for the initial filtered set, and must be
recomputed after mapping recovery before reporting FigS8 results. The selected
representative proteins, alignment definitions, and provenance are documented
above and in the retained mapping/alignment tables. The requested recovery
review does not authorize sequence-similarity partner selection, conservation
matching, adjustment of original outcomes, or changes to the fixed alignment
definitions.

## Measured universes and mapping exclusions

The recovery review completed on 2026-11-11. The user approved applying 628
evidence-supported mappings; 29 cases remain excluded. No mappings or new
alignments have yet been applied. The fresh mapping-only SLURM job is prepared
and awaits approval of its exact submission command. See
[`RECOVERY_REVIEW.md`](RECOVERY_REVIEW.md) for the per-gene evidence, residual
exclusions, and job specification.

| Category | Canonical eligible | Initially aligned/measured | One-to-one measured |
|---|---:|---:|---:|
| SSDOs | 568 | 527 | 101 |
| Disease excluding SSDOs | 4,512 | 4,195 | 750 |
| Orthologues excluding disease | 4,012 | 3,713 | 1,448 |

The initial filtered alignment has the same measured universe across all three
metrics: 8,435 genes
(2,299 one-to-one), with no additional alignment dropout. Report metric-specific
universes if downstream filtering changes this; do not silently force a shared
universe. All 8,435 retained genes had valid fly and human representatives.

Pre-alignment exclusions: 629 fly genes with multiple ENSG candidates, 18 with
unmapped HGNC IDs, nine SSDOs without supplied pairs, and one missing human
protein. The latter is FBgn0033337 / HGNC:13605 / ENSG00000125122. Separately,
906 LRG xref rows were audited as a different reference namespace, not counted
as ambiguous ENSG candidates. Of the 629 ambiguous observations, 618 had only one candidate with a valid
protein, ten had two, and one had none. These exclusions are under
reference-aware review; a candidate's protein availability alone is not
sufficient evidence to resolve a biological mapping.

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
- Exact direct HGNC Symbol xrefs map to ENSG IDs. LRG records are separately
  audited; zero/multiple ENSG candidates are exclusions. Validate the
  gene/transcript/protein chain against GTF and FASTA.
- For each species, choose a uniquely mapped valid Ensembl-canonical protein;
  otherwise choose the longest valid protein, breaking ties by stable protein
  ID. One fly and one human protein in the retained set used this fallback.
  Strip one terminal stop, retain versioned IDs, and exclude internal stops,
  conflicting mappings, and duplicate FASTA stable IDs.
- Exact global Gotoh alignment: Biopython 1.86, NumPy 1.26.4, BLOSUM62, affine
  gap cost `10 + (gap_length - 1) * 0.5`, including terminal gaps. Fly is the
  target; take the first optimal traceback from the pinned implementation,
  without choosing among ties by conservation. The exact matrix and parameters
  are in alignment provenance.
- Mask positions outside the 20 standard amino acids to X **for scoring only**.
  Original residues remain in saved alignments; no masked position contributes
  an identity match. Masking affects 297 pairs: 33 SSDOs, 155 disease-only,
  109 orthologue-only. It does not exclude these pairs.
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

| Job | Outcome |
|---|---|
| 4935248 preparation | Failed on an LRG namespace assumption, after all five reference files passed checksums; original scripts/logs retained |
| 4935338 preparation recovery | Completed, exit 0, 77 seconds; reused/reverified references; fresh mapping outputs |
| 4935733 global alignment | Completed, exit 0, 3m57s, about 1.36 GiB peak RSS; all 8,435 pairs aligned |

Historical CSVs, JSONs, logs, completion markers, and manifests were moved
**without rewriting bytes**. They intentionally record the original
`side/fly_human_conservation/` paths. Resolve that prefix to
`fly_human_conservation/` for the relocated files; `Data/` paths did not change.
Relative alignment paths inside the gene table are relative to
`alignment_4935733/`, and point to ignored HPC-only files.
The directory's `.gitattributes` preserves original CSV/TSV line endings,
empty fields, and archived provenance bytes during checkout.

The original executed source versions are in `provenance/executed_scripts/`;
their hashes match the successful jobs' provenance. The failed version is
in `failed_job_4935248_scripts/`. These are **archives, not runnable entry
points**. `provenance/HPC_EXECUTION_NOTES.md` preserves the full chronological
pre-move notes, with their original paths. Current entry scripts were updated
only for the new layout. No scientific calculation was rerun for the handoff.

## Tracked versus HPC-only

Git includes scripts/tests/requirements, compact tables, reference and output
manifests, completion records, original source snapshots, and successful/failed
logs. It excludes:

- `fly_human_conservation/references/`: compressed reference downloads.
- `fly_human_conservation/alignment_*/alignments/`: individual sequence/traceback
  JSONs (listed by checksum in the tracked output manifest).
- `fly_human_conservation/mapping_*/tables/protein_candidates.csv`: large
  protein-candidate audit, not needed for statistics on the approved gene set.
- `side/fly_human_conservation/.venv_alignment_biopython186_numpy1264/`:
  installed HPC environment, deliberately left at its original path because
  virtual environments are not portable.
- The original `Data/inputs/all_orthologs.csv` and
  `Data/inputs/all_disease_orthologs.csv`, still ignored. The tracked full-pair
  table already preserves their pair content and source membership.

After the compact refreshed tables are committed/pulled, the local R analysis
needs no untracked files. A raw-input rerun or an individual alignment
investigation requires the HPC-only references and alignment files. Never copy
the installed virtual environment to another machine; recreate it from
`alignment_requirements.txt` if needed.

## HPC reruns (not part of the local handoff)

Do not overwrite existing outputs. First inspect inputs/modules/logs, obtain
approval of the exact new submission, and create the log directory before
`sbatch`. Updated entry points, from the repository root, are:

```bash
mkdir -p fly_human_conservation/logs
sbatch fly_human_conservation/prepare_references.sh
sbatch fly_human_conservation/prepare_references.sh --reuse-verified-references
sbatch fly_human_conservation/align_proteins.sh
```

These are alternatives/stages, not a chain to execute now. Default preparation
refuses existing reference/table/provenance directories. Explicit reuse mode
checks unchanged raw input hashes and revalidates existing references, writing
a fresh `mapping_NEWJOBID/` directory. Alignment creates
`alignment_NEWJOBID/`; it checks all selected sequence hashes. For a refreshed
mapping, configure the alignment step to consume that specific new
`mapping_NEWJOBID/` output rather than the historical `mapping_4935338/`
directory. SLURM headers retain resources, email notifications, and new log
paths. The alignment wrapper uses the original ignored HPC environment path,
creating it from pinned requirements only if absent. Bulk references and raw
orthology inputs are needed for these HPC reruns, not for local summaries.
