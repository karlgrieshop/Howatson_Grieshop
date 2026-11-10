# Fly-human conservation: local analysis handoff

**HPC preparation and global alignment are complete.** This remains exploratory,
separate from canonical gene sets and manuscript figures/tables. On 2026-11-10,
the user requested moving the work from ignored `side/fly_human_conservation/`
to this tracked directory. No sequences, metrics, categories, or exclusion rules
were changed during that move.

**Everything needed for local statistics and plotting comes through Git.
No `scp`, ignored raw inputs, reference downloads, alignment files, or Python
environment are required for those next steps.** The user handles push/pull.

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

## Local agent: next analysis only

1. Read `.github/copilot-instructions.md` and this README. Do not rerun
   preparation/alignment or change canonical analyses.
2. Read the main gene table directly. Compare the three mutually exclusive
   categories in the order SSDOs; disease excluding SSDOs; orthologues excluding
   disease. Do not add Background or impute excluded genes as zero.
3. For each of `identity_percent`, `fly_coverage_percent`, and
   `human_coverage_percent`, export group arithmetic means and 95% percentile
   bootstrap intervals (5,000 resamples, seed 1). Plot the observed arithmetic
   means with intervals, measured sample sizes, and explicit percentage axes.
4. For each metric separately, run the three two-sided pairwise differences
   in arithmetic means (10,000 two-sample permutations; reproducible documented
   seeds). Use the existing repository convention
   `(exceedances + 1) / (n_perm + 1)` and apply BH across the three contrasts
   within that metric. Keep the two directional coverage metrics distinct.
5. Repeat on the existing one-to-one subset as a sensitivity analysis, with
   its own reported denominators and separately labelled per-metric BH families.
   Export raw/adjusted P-values, group sizes, mean differences, seeds, counts,
   method, and bootstrap/plot-data tables. Do not interpret non-significance as
   equivalence; shared human partners also limit independence in the main set.
6. Follow the established R style. The relevant prior art is
   `Plotting_and_analyses/boot_permute_metrics.R`; do not source its entire
   analysis driver, which requires unrelated ignored inputs and writes canonical
   results. Keep new scripts, figures, and summaries inside
   `fly_human_conservation/`, e.g. `summaries/` and `figures/`.

Statistical tests, confidence intervals, and plots have **not** been run yet.
The user explicitly chose to leave potentially resolvable mapping exclusions
excluded. Do not rescue them, match on conservation, adjust original outcomes,
or change representative proteins in this follow-up.

## Measured universes and preliminary means

| Category | Canonical eligible | Aligned/measured | One-to-one measured |
|---|---:|---:|---:|
| SSDOs | 568 | 527 | 101 |
| Disease excluding SSDOs | 4,512 | 4,195 | 750 |
| Orthologues excluding disease | 4,012 | 3,713 | 1,448 |

All three metrics currently have the same measured universe: 8,435 genes
(2,299 one-to-one), with no additional alignment dropout. Report metric-specific
universes if downstream filtering changes this; do not silently force a shared
universe. All 9,092 eligible fly genes had valid representatives.

Pre-alignment exclusions: 629 fly genes with multiple ENSG candidates, 18 with
unmapped HGNC IDs, nine SSDOs without supplied pairs, and one missing human
protein. The latter is FBgn0033337 / HGNC:13605 / ENSG00000125122. Separately,
906 LRG xref rows were audited as a different reference namespace, not counted
as ambiguous ENSG candidates. Of the 629 ambiguous observations, 618 had only
one candidate with a valid protein, ten had two, and one had none; all remain
excluded by the user's decision.

| Category | Mean identity (%) | Mean fly coverage (%) | Mean human coverage (%) |
|---|---:|---:|---:|
| SSDOs | 34.9472 | 83.5753 | 88.7848 |
| Disease excluding SSDOs | 35.7591 | 85.6610 | 89.1863 |
| Orthologues excluding disease | 34.9032 | 84.9744 | 88.5492 |

These are descriptive arithmetic means, not significance/equivalence results.

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
Fly: BDGP6.54 / Ensembl release 115, using the peptide and GTF files named in
`HPC_AGENT_PROMPT.md`. Original URLs, retrieval dates, provider MD5/BSD checks,
and SHA-256 fingerprints are retained in reference manifests.

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

**No untracked file is required for the local next analysis.** Only a raw-input
rerun or an individual sequence-alignment investigation would need transfers or
reference retrieval. Never copy the installed virtual environment to another
machine; recreate it from `alignment_requirements.txt` if needed.

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
`alignment_NEWJOBID/`; it checks all selected sequence hashes and uses the
fixed successful `mapping_4935338/` inputs. SLURM headers retain resources,
email notifications, and new log paths. The current alignment wrapper uses the
original ignored HPC environment path, creating it from pinned requirements
only if absent. Bulk references and raw orthology inputs are needed for these
HPC reruns, not for local summaries.
