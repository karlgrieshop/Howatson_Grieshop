# Fly-human protein conservation: exploratory preparation

Submit from the repository root, after approval of this exact command:

```bash
mkdir -p side/fly_human_conservation/logs
sbatch side/fly_human_conservation/prepare_references.sh
```

This stage retrieves fixed-release references and audits gene/protein mappings.
It does not align sequences, compute conservation, run tests, or alter canonical
analyses. Python 3.12.7 (`python/anaconda/2024.10/3.12.7`) uses only its standard
library and the existing GNU `sum` utility. Resources: compute partition, bio
account, one CPU, 8 GB RAM, two hours. SLURM email matches the supplied examples.

Preflight on 2026-11-10 found 568 SSDOs, 4,512 disease-only genes, and 4,012
orthologue-only genes. Supplied pairs cover 559, 4,512, and 4,012 respectively;
the nine missing SSDOs remain in the audit. Full input: 35,107 unique pairs,
35,047 within the canonical universe; all 10,574 disease-input pairs are a
score-consistent subset. Scores span 3-14. Full-graph one-to-one candidates
number 104, 796, and 1,552 before reference mapping. These are preflight counts,
not conservation sample sizes. The approved preparation tolerates explicitly
reported gene dropout without changing categories or rescuing lower-score pairs.

An exact global aligner has not been selected: checked Python environments lack
Biopython/parasail, and R environments lack Biostrings/pwalign. MAFFT 7.526 is
available but is not silently substituted. No software is installed by this job.

Approved submission on 2026-11-10: job **4935248**, initially **PENDING** in
`compute`, using the command above. Acceptance is not evidence of success.
No ongoing monitoring or dependent submission was started.

## Failure diagnosis and proposed recovery

Job 4935248 failed after 73 seconds with exit code 1, using about 85 MiB peak
resident memory. All five bulk reference files passed provider checksum checks;
the failure occurred during xref parsing, not download or resource allocation.
The parser incorrectly treated the valid `LRG_10` record for `HGNC:9255` (PPIB)
as a malformed ENSG gene ID. LRG records are a separate reference namespace,
not an additional ENSG candidate and not evidence that the HGNC gene is absent.

The corrected parser audits LRG rows in `excluded_human_xrefs.csv` and uses only
ENSG records for the Ensembl-gene mapping. Unknown namespaces still fail with a
diagnostic. No category, DIOPT ranking, or representative-protein rule changed.
The original scripts and README are archived in `failed_job_4935248_scripts/`;
their entry-script hashes match the failed job's `provenance/execution.json`.
Original references, provenance, empty tables directory, and failed logs remain
untouched.

Recovery command explicitly approved and submitted on 2026-11-10:

```bash
sbatch side/fly_human_conservation/prepare_references.sh --reuse-verified-references
```

Recovery job **4935338** completed with exit code 0 in 77 seconds, using about
557 MiB peak resident memory. Outputs: `mapping_4935338/{tables,provenance}/`;
logs: `logs/prepare-4935338.{out,err}`. Stderr is empty and the completion marker
exists. Completion was inspected on 2026-11-10 after the user reported it;
there was no continuous monitoring. Retain the successful and failed logs.

Resources, modules, email notifications, and log paths are unchanged. This mode
reads the original `references/`, `provenance/reference_manifest.tsv`, and
`provenance/input_manifest.csv`. It requires unchanged input fingerprints and
revalidates every reference's size/SHA-256 plus pinned provider checksums for
the five bulk files. It does not download or overwrite original results.

Outputs go into a new `mapping_NEWJOBID/{tables,provenance}/` directory.
Copied reference manifests retain the original retrieval dates and paths
relative to the source directory recorded in `execution.json`, adding a
`reuse_verified_utc` timestamp. Recovery refuses any existing destination.
Successful preparation still requires its own `PREPARATION_COMPLETE.json`.

## Verified preparation results

| Category | Canonical eligible | Selected pairs | Mapping-ready | One-to-one mapping-ready |
|---|---:|---:|---:|---:|
| SSDOs | 568 | 559 | 527 | 101 |
| Disease excluding SSDOs | 4,512 | 4,512 | 4,195 | 750 |
| Orthologues excluding disease | 4,012 | 4,012 | 3,713 | 1,448 |

These are preparation counts, not measured identity/coverage universes.
There are 8,435 mapping-ready fly genes of 9,092 eligible (2,299 one-to-one);
they use 6,222 distinct human proteins. Exclusions: 629 fly genes with multiple
ENSG candidates, 18 with unmapped HGNC IDs, nine without supplied pairs, and
one without a valid human protein. All 9,092 fly genes have valid representatives.
The missing human-protein pair is FBgn0033337 / HGNC:13605 / ENSG00000125122.
LRG exclusions comprise 906 xref rows, separately audited rather than counted
as ambiguous ENSG mappings.

Among the 629 ambiguous fly-gene observations, 618 have exactly one candidate
with a valid protein in these references, ten have two, and one has none.
This may permit reference-aware resolution after reviewing locus/annotation
status, without sequence-similarity selection; no mapping rule was changed and
these observations remain excluded. One fly and one human representative among
the mapping-ready set use the longest-valid fallback; the others use canonical
proteins. Nonstandard residues occur in 287 unique fly representatives and six
unique human representatives, requiring an explicit alignment/scoring policy.

Gene IDs are unique in the 9,092-row audit and 9,083-row selected table;
mapping-ready protein joins and category denominators agree with the summary.
At preparation completion, no conservation metrics, bootstrap summaries,
permutation tests, or plots had been computed.

## Approved alignment stage

The user directed that ambiguous/excluded mappings stay excluded for now.
Proceed only with the 8,435 mapping-ready pairs from `mapping_4935338/`.
Compact-table inspection found nonstandard residues in 297 pairs: 33 SSDOs,
155 disease-only, and 109 orthologue-only. The largest sequence-length product
is 18,468 fly residues by 35,991 human residues; total products across retained
pairs sum to 6,461,437,817 dynamic-programming cells. No alignments were run to
obtain these counts.

Approved on 2026-11-10 and implemented in `align_proteins.py`: exact global alignment with
Biopython PairwiseAligner, BLOSUM62, affine gap opening 10/extension 0.5,
penalized terminal gaps (length-k gap costs 10 + (k-1)*0.5). The pinned library
uses the Gotoh global alignment algorithm. With fly as target and human as
query, take the first optimal traceback returned by Biopython; do not enumerate
all tied optima or select among them using conservation. Save the exact matrix,
all terminal/internal gap scores, score tolerance, and package versions.
Source: <https://biopython.org/docs/1.86/Tutorial/chapter_pairwise.html>.

Installed approved Biopython 1.86 and NumPy 1.26.4 from binary wheels in
`.venv_alignment_biopython186_numpy1264/` using Python 3.12.7 from
`python/anaconda/2024.10/3.12.7`. System site packages are disabled; no shared
software was modified. Pins are in `alignment_requirements.txt`; the wheel
URLs/hashes are in the environment's `install_report.json`. Environment setup
and toy validation were lightweight; all project alignments run through SLURM.
The entry script creates this isolated environment if absent, or checks the
existing dependency versions and runs toy tests before biological computation.
The user approved the nonstandard-residue masking policy on 2026-11-10:
mask positions outside the 20 standard amino acids to X for scoring only,
retain original sequence positions, and never count masked positions as exact
identity matches. This preserves the 297 affected pairs rather than dropping
them. Per-protein residue counts and all original sequences remain auditable.

Metrics: exact standard-residue matches divided by all alignment
columns, and separate fly/human coverage as residue-residue columns divided
by the respective full protein length, each reported in percent. Preserve
original sequences, aligned sequences, scores, masking counts, sequence hashes,
and per-gene numerators/denominators. Use one fixed scoring scheme for all pairs;
never select orthologues or isoforms by alignment results.
Masked positions count as residues for coverage, but never as identity matches.
These are whole-alignment identity and paired-residue coverage, not BLAST HSP
identity or percent span covered (which would be 100% in a global alignment).

Submit the approved command from the repository root:

```bash
sbatch side/fly_human_conservation/align_proteins.sh
```

Submitted on 2026-11-10 as job **4935733**, initially **PENDING** in `compute`.
After the user reported completion, accounting confirmed **COMPLETED**, exit
0, elapsed 3m57s, peak resident memory about 1.36 GiB. Outputs:
`alignment_4935733/`; logs: `logs/align-4935733.{out,err}`. Stderr contains
successful toy-test output, not an alignment error. No continuous monitoring
or dependent submission was performed; all successful/failed logs remain.

Resources: `compute`, account `bio`,
one CPU, 16 GB RAM, eight hours; same email notifications as preparation.
Runtime remains uncertain without a compute-node benchmark. Inputs: the
successful preparation's selected-pair/representative tables and completion
provenance, plus the original checksum-verified fly/human peptide references.
Outputs: new `alignment_JOBID/{alignments,tables,provenance}/`, with logs
`logs/align-JOBID.{out,err}`. No overwriting or reference changes. This stage
produces metrics only; bootstrap/permutation summaries and plotting are
a later separately reviewed step.

- `tables/gene_conservation.csv`: one measured row per retained fly gene,
  category, original pair, IDs, selection rules, lengths, exact numerators and
  denominators, scores, three percentages, masking counts, hashes, alignment path.
- `tables/one_to_one_conservation.csv`: subset using the existing full-graph
  one-to-one flag, not a new post-alignment definition.
- `tables/category_counts.csv`: eligible, pre-alignment excluded, aligned,
  per-metric measured, masked, and one-to-one counts by category. Positive-length
  global alignments yield all three metrics for the same retained genes here;
  separate count columns disclose this rather than forcing missing metrics.
- `alignments/FBgnID.json`: raw and masked sequences, original residues in
  aligned coordinates, optimal traceback coordinates, IDs, and score.
- `provenance/`: input/script hashes, versions, scoring matrix/parameters,
  output hashes, and `ALIGNMENT_COMPLETE.json` written last. Missing completion
  means partial results, not success. Unexpected data/alignment errors stop
  explicitly, rather than silently dropping genes.

Preflight checks canonical/orthology fingerprints against successful preparation,
selected/audit/representative consistency, measured-universe counts against its
completion marker, and peptide reference hashes against its verified manifest.
Each selected sequence must match its original ID, transcript, length, terminal
stop treatment, nonstandard count, and sequence hash. Source files are rechecked
before the completion marker, so edits during execution fail explicitly.
Only small synthetic examples are used for login-node tests. Keep all logs,
partial outputs, and original preparation directories.

### Verified alignment results (2026-11-10)

All 8,435 retained pairs aligned; there was no additional alignment dropout.
All three metrics have the same measured universe here: 527 SSDOs / 568
canonical eligible, 4,195 disease-only / 4,512 eligible, and 3,713 orthologue-only
/ 4,012 eligible. The one-to-one subset contains 101, 750, and 1,448 genes
respectively (2,299 total); 297 aligned pairs have masked positions.

| Category | Measured n | Mean identity (%) | Mean fly coverage (%) | Mean human coverage (%) |
|---|---:|---:|---:|---:|
| SSDOs | 527 | 34.9472 | 83.5753 | 88.7848 |
| Disease excluding SSDOs | 4,195 | 35.7591 | 85.6610 | 89.1863 |
| Orthologues excluding disease | 3,713 | 34.9032 | 84.9744 | 88.5492 |

These are descriptive arithmetic means only, not significance/equivalence
claims. Masked residues remain in the full-length denominators but contribute
no identity matches. Bootstrap intervals, permutation tests, BH corrections,
one-to-one statistical sensitivity, and plots have not been run.

Verified completion-manifest/table hashes, unique gene IDs and sensitivity
membership, metric numerators/denominators and finite 0-100 ranges for every
gene, and one alignment file per measured gene. Independently reconstructed
metrics for three saved alignments, including the largest and a masked pair.
The next analysis stage is group bootstrap/permutation summaries and figures
using these compact CSVs; no alignment rerun or bulk-reference transfer is
needed for that stage.

## Inputs and references

- `Data/inputs/all_orthologs.csv`, `Data/inputs/all_disease_orthologs.csv`.
- `Data/outputs/sexspecific_disease_ids.csv`,
  `Data/outputs/all_disease_ortholog_ids.csv`,
  `Data/outputs/all_ortholog_ids.csv`, `Data/outputs/dmel_gene_reference.csv`.
- Human GRCh38.p14, GCA_000001405.29, genebuild 2026-04:
  <https://ftp.ebi.ac.uk/pub/ensemblorganisms/GCA/000/001/405/29/ensembl/2026_04/geneset/>.
  Retrieve `pep.fa.bgz`, `genes.gtf.bgz`, `xref.tsv.gz`, and `md5sum.txt`.
- Fly BDGP6.54 / Ensembl release 115:
  <https://ftp.ensembl.org/pub/release-115/fasta/drosophila_melanogaster/pep/Drosophila_melanogaster.BDGP6.54.pep.all.fa.gz>
  and
  <https://ftp.ensembl.org/pub/release-115/gtf/drosophila_melanogaster/Drosophila_melanogaster.BDGP6.54.115.gtf.gz>,
  with each directory's `CHECKSUMS`.

The script pins the checksums observed on 2026-11-10, verifies both the provider
manifest and downloaded bytes, and records SHA-256 fingerprints. Human checks
use MD5; fly checks use BSD `sum -r` and 1-KiB block counts. BGZF files are read
as concatenated gzip members without modifying the downloads.

## Definitions and exclusions

Canonical categories are SSDOs, disease IDs minus SSDOs, and orthologue IDs
minus disease IDs. All must lie in the tracked protein-coding universe.
Background is excluded, never assigned zero conservation.

Retain the complete input pair table, including pairs outside that universe.
The disease table must be a score-consistent subset of the full table. For each
canonical fly gene, select its highest-score pair from the full table; break
ties lexicographically by stable HGNC ID. No new DIOPT threshold, similarity
selection, disease-partner restriction, or lower-score mapping rescue is used.
The observed input minimum of 3 does not establish upstream filtering provenance.

Map HGNC IDs only through exact `HGNC Symbol`, `DIRECT` ENSG xrefs. Audit LRG
reference records separately rather than counting them as ENSG candidates.
Zero or multiple distinct ENSG mappings are explicit exclusions. Strip stable-ID version suffixes
for joins, retaining the selected versioned FASTA protein ID. Parse the GTF
gene/transcript/protein chain and confirm it against FASTA headers.

For both species, select a unique valid Ensembl-canonical protein when available;
otherwise choose the longest valid protein, breaking ties by stable protein ID.
Multiple canonical candidates trigger the recorded longest-valid fallback.
Valid proteins are nonempty with no internal stops or unsupported characters.
Allow the 20 standard amino acids plus B, X, Z, J, U, and O, reporting the count
of these additional residues. Strip and record one terminal stop. Duplicate
FASTA stable IDs and conflicting gene/transcript mappings are excluded, not
silently resolved. Nonstandard residues require an explicit policy before
alignment; this preparation does not imply they are scoreable by an aligner.

One-to-one sensitivity requires degree one at both ends of the full supplied
bipartite graph, before canonical filtering, selection, or mapping exclusions.
Tables distinguish eligible, selected, mapped, and alignment-ready genes;
`ready_for_alignment` means mapping/sequence preparation only, not a measured
conservation result. Reasons may overlap; only total excluded genes is disjoint
from ready genes. No-pair genes stay in the canonical audit.

## Outputs and safeguards

Everything remains under this ignored side directory:

- `references/{human,fly}/`: original compressed references and checksums.
- `provenance/reference_manifest.tsv`: retrieval dates, URLs, provider checks,
  sizes, SHA-256 fingerprints (tab-delimited).
- `provenance/input_manifest.csv`, `execution.json`, and annotation headers:
  input hashes, script hashes, Git revision, Python version, job ID, and builds.
- `tables/all_pairs.csv`, `selected_pairs.csv`, `one_to_one_pairs.csv`:
  pair provenance, selection, and mapping status.
- `tables/hgnc_mapping_audit.csv`, `protein_candidates.csv`,
  `representative_proteins.csv`, `gene_mapping_audit.csv`: complete mapping
  candidates, protein choices, and gene-level exclusions.
- `tables/excluded_human_xrefs.csv`: LRG xref rows with HGNC IDs and reasons;
  the HGNC audit separately counts retained ENSG and excluded LRG rows.
- `tables/diopt_distributions.csv`, `category_attrition.csv`: score distributions
  and category-specific denominators and missingness.
- `logs/prepare-JOBID.{out,err}`: retain successful and failed logs.
- `provenance/PREPARATION_COMPLETE.json`: written only after all exports and
  unchanged-input checks succeed. Its absence means preparation is incomplete.

In default download mode, existing reference/table/provenance directories cause
an immediate failure. The explicitly approved reuse mode described above writes
only a fresh job-specific directory. Neither mode automatically overwrites,
cleans up, or retries. Partial outputs remain for diagnosis. Inspect only when the user requests completion/failure
follow-up; do not monitor a submitted job continuously.

After output verification, transfer the compact `tables/`, `provenance/`, this
README, and both entry scripts to the same ignored path in the local clone.
Bulk references need not leave the HPC. Provide a user-run `scp` command once
the external SSH hostname and local destination are known; never automate
authentication. Alignment method, residue handling, conservation definitions,
and any subsequent SLURM command require separate review/approval.
