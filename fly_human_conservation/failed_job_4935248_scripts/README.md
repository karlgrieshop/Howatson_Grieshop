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

Map HGNC IDs only through exact `HGNC Symbol`, `DIRECT` xrefs. Zero or multiple
distinct ENSG mappings are explicit exclusions. Strip stable-ID version suffixes
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
- `tables/diopt_distributions.csv`, `category_attrition.csv`: score distributions
  and category-specific denominators and missingness.
- `logs/prepare-JOBID.{out,err}`: retain successful and failed logs.
- `provenance/PREPARATION_COMPLETE.json`: written only after all exports and
  unchanged-input checks succeed. Its absence means preparation is incomplete.

Existing reference/table/provenance directories cause an immediate failure:
there is no automatic overwrite, cleanup, retry, or resume. Partial outputs
remain for diagnosis. Inspect only when the user requests completion/failure
follow-up; do not monitor a submitted job continuously.

After output verification, transfer the compact `tables/`, `provenance/`, this
README, and both entry scripts to the same ignored path in the local clone.
Bulk references need not leave the HPC. Provide a user-run `scp` command once
the external SSH hostname and local destination are known; never automate
authentication. Alignment method, residue handling, conservation definitions,
and any subsequent SLURM command require separate review/approval.
