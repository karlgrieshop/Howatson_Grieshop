# Mapping recovery review: 2026-11-11

This is a review of the original 657 exclusions, not a change to canonical
categories, original pair ranking, protein rules, alignment definitions, or
FigS8 statistics. Original mapping/alignment directories and logs are untouched.
No recovery job has been submitted.

## Initial per-gene review

Read `recovery_review_preflight_20261111/tables/per_gene_resolution.csv`:
657 unique fly-gene rows, with original reason, selected human partner/symbol
and score, candidate ENSGs, original one-to-one flag, original disease-source
records, evidence/source/release, decision, proposed IDs, and uncertainty.

This is explicitly **preflight**, not a completed full-reference adjudication.
Only compact tracked audits and two bounded HGNC REST records were read.
The HGNC JSONs, retrieval URLs/dates, hashes, input/output manifests, and
`REVIEW_COMPLETE.json` are under the preflight directory's `provenance/`.
No original exclusion is yet recovered/applied. Category counts remain:

| Category | Canonical eligible | Originally measured | Excluded | Full-graph one-to-one eligible / measured / excluded |
|---|---:|---:|---:|---|
| SSDOs | 568 | 527 | 41 | 104 / 101 / 3 |
| Disease excluding SSDOs | 4,512 | 4,195 | 317 | 796 / 750 / 46 |
| Orthologues excluding disease | 4,012 | 3,713 | 299 | 1,552 / 1,448 / 104 |

`tables/category_attrition_review.csv` separates original measured counts,
proposed additions, projected counts if approved, and residual exclusions.
Proposed counts are not measured counts; no fresh alignments have run.

## Findings and recovery criteria

### Multiple ENSG candidates: 629 fly genes, 388 distinct HGNC partners

Current HGNC explicitly assigns HGNC:6698 (LRP6) to ENSG00000070018, one of
the candidates for FBgn0000119. This is an authoritative identity lead, not
recovery based only on protein availability. The fixed-reference gene name,
chromosome, coding status, direct HGNC xref, and representative protein still
need corroboration in the full review.

The proposed job retrieves the official current HGNC complete-set TSV,
linked from <https://www.genenames.org/download/archive/>. It records the
snapshot date/hash; current identity evidence is distinct from the fixed
Ensembl 2026-04 sequence/annotation release, which does not change.

For a proposed ambiguous resolution, require all of:

- The unchanged selected HGNC ID is Approved and designated protein-coding,
  with one explicit ENSG ID in its authoritative HGNC record.
- That ENSG is an original candidate, lies on a main reference chromosome,
  and exists in the fixed GTF as protein-coding.
- Its fixed GTF name matches the HGNC approved or previous approved symbol;
  its chromosome agrees with HGNC's cytogenetic location.
- Fixed-reference DIRECT HGNC Symbol xrefs support that HGNC ID without a
  competing direct HGNC assignment for the selected ENSG.
- Only after identity is supported, the unchanged representative-protein
  rules yield a valid protein in the same fixed peptide release.

Export the annotation and direct xrefs for **every** original candidate,
including alternatives absent from the fixed GTF. A protein's availability,
length, identity, or alignment score must never decide biological identity.
Any discrepancy or incomplete evidence stays manual/unresolved; no proposal
is automatically applied.

### Unmapped HGNC: 18 fly genes

All 18 source IDs are **HGNC:0**, a placeholder shared across multiple human
symbols, not 18 valid withdrawn HGNC records. Their source symbols are:
`LOC102724770`, `CRYAA2`, `U2AF1L5`, `CBSL`, `LOC400499`,
`LOC102724159`, `TBC1D7-LOC100130357`, and `LOC102724488`.

The full review checks exact currently approved HGNC symbol assignments,
preserving each selected source partner and score. Only an exact unique
approved symbol plus the reference corroboration above can produce a proposal.
Previous/alias symbol matches and fixed-reference NCBI Gene xrefs (including
the numeric ID in a LOC label) are exported as history leads, not silently
equated to another gene. Unresolved cases may need a targeted NCBI/history
follow-up after examining these records.

The original full-graph one-to-one flags remain unchanged, including their
dependence on the original placeholder namespace. The review does not
recompute graph degrees using repaired human identifiers. Do not claim that
the original HGNC:0 grouping represents one biological gene.

### No valid protein: FBgn0033337

HGNC:13605 explicitly maps to ENSG00000125122 but is now **FBXL9P, a pseudogene**
(HGNC record modified 2026-04-30; previous symbols FBXL9/LRRC29). This supports
retaining exclusion from a protein-sequence comparison unless fixed-release
annotation provides valid, consistent contrary evidence warranting review.
Do not use a historical protein, synthesize a translation, change release, or
substitute another human orthologue. The preflight marks it
`retain_excluded_HGNC_nonprotein_locus`, not a recoverable missing download.

Source: <https://rest.genenames.org/fetch/hgnc_id/HGNC:13605>.

### No supplied pairs: nine SSDOs

The tracked `Data/outputs/sexspecific_disease_gene_data.csv` identifies:

| Fly gene | Symbol | Retained source evidence |
|---|---|---|
| FBgn0003391 | shg | FlyBase prostate-cancer search |
| FBgn0015609 | CadN | FlyBase prostate-cancer search |
| FBgn0035670 | CG10472 | Literature source `alzheimer_disease_orthologs.csv`, human symbol CTRC |
| FBgn0036287 | CG10663 | FlyBase hemophilia search |
| FBgn0036612 | CG4998 | Literature source `coronary_artery_disease_orthologs.csv`, human symbol F12 |
| FBgn0050169 | Brca2 | FlyBase breast-cancer, gonadal/ovarian-dysgenesis, and prostate-cancer searches |
| FBgn0052772 | CG32772 | FlyBase retinitis-pigmentosa search |
| FBgn0260450 | CalpC | FlyBase type-2-diabetes search |
| FBgn0261679 | CG42726 | FlyBase retinitis-pigmentosa search |

`Data/inputs/ssdo_inputs/` is absent on this HPC clone. The compiler retains
the CTRC/F12 symbol provenance, but not the original pair scores/ranks.
For literal original records, recover these two source CSVs from
`Data/inputs/ssdo_inputs/orthology_sources/`, plus the seven distinct text
files named in the table under `Data/inputs/ssdo_inputs/sex_specific_diseases/`:
`prostate_cancer.txt`, `hemophilia.txt`, `breast_cancer.txt`,
`gonadal_dysgenesis.txt`, `ovarian_dysgenesis.txt`,
`retinitis_pigmentosa.txt`, `type2_diabetes.txt`.
The original FlyBase query results may contain only fly IDs, so authoritative
FlyBase disease-model/orthology records may still be needed to establish
human partners.

Do not guess BRCA2/CDH partners from fly symbol resemblance, infer orthologues
from protein similarity, invent a DIOPT score, or claim uniqueness from one
retained human-symbol row. Any proposed added pair requires separate source
review and approval. The review leaves these unpaired, with one-to-one status
explicitly unassigned rather than treating the audit's False as proof of
many-to-many orthology.

## Proposed review-only SLURM job: awaiting approval

```bash
sbatch fly_human_conservation/review_mapping_recovery.sh
```

Resources in headers: `compute`, account `bio`, one node/task/CPU, **8 GB**,
**2 hours**, job `fly_human_review`; email `ALL` to the existing address.
Log directory already exists; output/error:
`fly_human_conservation/logs/review-%j.{out,err}`.
Module: verified `python/anaconda/2024.10/3.12.7`; standard library only.
The partition/account and original successful logs were inspected.

Required inputs:

- Original `mapping_4935338/tables/` gene/pair/HGNC/xref/category audits, and
  original preparation/alignment manifests/completion records used by
  `check_handoff.py`.
- `Data/outputs/sexspecific_disease_gene_data.csv` and the tracked canonical
  ID/reference tables. The current pair CSVs under `Data/inputs/` are present;
  existing selected pairs are authoritative and are not re-ranked.
- All eight downloaded reference/manifest files under
  `fly_human_conservation/references/{human,fly}/`, as listed by the original
  reference manifest. Reverify existing fingerprints/provider checksums.
- Network access to
  `https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt`
  (about 17 MB; HEAD check succeeded). Snapshot this official current HGNC
  identity dataset only; do not download another sequence/annotation release.

New output directory: `fly_human_conservation/recovery_review_NEWJOBID/`.
It refuses an existing destination, and writes:

- `tables/per_gene_resolution.csv`: all 657 exclusions, with evidence,
  proposed IDs/decisions, and residual uncertainty.
- `tables/category_attrition_review.csv`: original and projected—not applied—
  category/one-to-one counts.
- `tables/candidate_gene_evidence.csv`, `candidate_reference_xrefs.csv`,
  `authoritative_HGNC_records.csv`: compact source identity/annotation audit.
- `tables/candidate_proteins.csv`: ignored per-protein validity audit.
- `downloads/hgnc_complete_set.tsv`: ignored immutable downloaded snapshot.
- `provenance/`: input/reference/source hashes, dates, source documentation,
  script/runtime/Git identity, output manifest, and `REVIEW_COMPLETE.json`.

Network, schema, checksum, annotation, and input-change errors stop explicitly;
partial files/logs are retained. The successful baseline directories, local
FigS8 summaries, and all old logs remain unchanged. No mapping is applied and
no alignments or statistical analyses run in this review job.

After review, report supported recoveries and residual exclusions, and obtain
approval for a separate fresh `mapping_NEWJOBID/` and `alignment_NEWJOBID/`
generation using unchanged selection/protein/alignment definitions. Preserve
original flags; nine genuinely added pairs need an approved graph-sensitivity
policy rather than silently redefining the existing one-to-one universe.
Track compact tables/scripts/provenance only; do not push.
