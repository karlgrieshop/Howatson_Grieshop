# Mapping recovery review: 2026-11-11

This is a review of the original 657 exclusions, not a change to canonical
categories, original pair ranking, protein rules, alignment definitions, or
FigS8 statistics. Original mapping/alignment directories and logs are untouched.
The review job has completed, but no mapping changes have been applied and no
new alignments have been run.

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

## Completed review-only SLURM job: 4955170

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

Job `4955170` (`fly_human_review`) completed with exit code 0 in 1:23, using
230,304 KiB maximum RSS. Its output manifest verifies all 11 listed files.
The fixed human GRCh38.p14/Ensembl 2026-04 and fly Ensembl 115 reference
fingerprints were revalidated; the current HGNC complete-set TSV was recorded
as a separate, newer identity snapshot. The job did not apply mappings or run
alignments.

### Results

The per-gene results are in
`recovery_review_4955170/tables/per_gene_resolution.csv` (657 unique genes).
The supporting compact audits are `candidate_gene_evidence.csv`,
`candidate_reference_xrefs.csv`, and `authoritative_HGNC_records.csv`; all are
alongside the table under `tables/`. Hashes, inputs, source snapshot, reference
checks, execution, and completion record are under `provenance/`.
The 17 MB downloaded HGNC snapshot and per-protein audit remain ignored bulk
artifacts; they are hashed in the output manifest.

The review proposes 628 mappings: 31 SSDOs, 303 disease-only orthologues, and
294 orthologue-only genes. Each proposal is grounded in the selected,
unchanged HGNC partner's Approved protein-coding record and unique ENSG
assignment; that ENSG is an original candidate, and the fixed annotation
corroborates coding status, current/previous approved symbol, chromosome, and
unique direct HGNC xref. A valid representative protein was then confirmed
under the original rule. This is stronger evidence than protein availability
alone, but current HGNC identity records postdate the fixed annotation
genebuild. The user approved applying exactly these 628 proposals on
2026-11-11; they have not yet been applied.

The 29 residual exclusions are:

- 18 `HGNC:0` placeholder cases: the current snapshot did not yield an exact,
  unique approved-symbol assignment passing the reference checks. The
  placeholder cannot establish partner identity; retain these as unresolved.
- Nine SSDOs without supplied pairs: original source files are absent, and the
  retained CTRC/F12 symbol leads do not preserve partner ranking or DIOPT
  scores. No pairs are invented.
- Two non-protein-coding loci: FBgn0033337 → HGNC:13605 / ENSG00000125122
  (FBXL9P pseudogene), and FBgn0036997 → HGNC:15547 / ENSG00000154316 (TDH
  locus classified non-protein-coding). Neither is substituted with a
  historical or alternative protein.

Projected mapped counts after the approved recoveries are 558/568 SSDOs,
4,498/4,512 disease-only, and 4,007/4,012 orthologue-only (9,063/9,092 total).
Among these proposals, 3, 46, and 103 respectively carry the original
full-graph one-to-one flag; these flags have not been recomputed. The
projections are not measured counts; no fresh alignments have run.

## Approved mapping job: awaiting exact-command approval

The mapping-only runner is `fly_human_conservation/recover_mapping.sh`. Before
submission it verifies the exact approved review-table hash, all 657 excluded
genes, all 628 selected-pair/candidate/protein correspondences, unchanged
category and original one-to-one flag, and the expected projected attrition.
Only the 628 reviewed decisions marked
`propose_recovery_authoritative_identity` are applied. The remaining 29 are
retained as excluded; no pair selection, category, reference, protein, or
alignment definition changes.

Proposed command, **not yet submitted**:

```bash
sbatch fly_human_conservation/recover_mapping.sh
```

- Resources: `compute`, account `bio`, 1 CPU, 8 GB, 2 hours; email notifications
  to the existing address.
- Inputs: canonical/source inputs and verified reference files used by
  `mapping_4935338`, plus the tracked approved
  `recovery_review_4955170/tables/per_gene_resolution.csv`. No external
  downloads.
- Outputs: fresh `mapping_JOBID/` with updated mapping tables, per-gene
  recovery audit, input/reference/output manifests, and completion record;
  logs `logs/recover-map-JOBID.{out,err}`. It refuses to overwrite an existing
  job directory.
- This job does **not** align proteins. After mapping completes and is
  inspected, a separate alignment command and resource request will be
  presented for approval.

The exact command still requires user approval before submission. Track compact
tables/scripts/provenance only; do not push.
