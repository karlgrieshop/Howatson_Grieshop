# HPC handoff prompt: fly–human protein conservation

Work in the existing repository clone at:

`/gpfs/data/GrieshopLab/Karl/SSDOs_Howatson_Grieshop`

Follow `.github/copilot-instructions.md`, especially its scientific safeguards
and SLURM submission protocol. This is a new, exploratory analysis. Keep its
scripts, logs, intermediate files, and outputs under the ignored
`side/fly_human_conservation/` directory. Do not change canonical gene sets,
existing figures/tables, or manuscript analyses.

## Objective

Assess whether average fly–human protein sequence identity and alignment
coverage differ among these existing, mutually exclusive categories:

1. SSDOs
2. Disease orthologues excluding SSDOs
3. Orthologues excluding all disease-associated genes

This is the initial conservation comparison only. Do not attempt conservation
matching or covariate adjustment of the original outcomes in this task.
Background genes have no retained human orthologue, so do not assign them zero
identity/coverage or include them in these direct orthologue comparisons.

## Start-up and inputs

1. From the repository root, inspect `git status`, relevant `Data/inputs` and
   `Data/outputs` files, available modules/tools, SLURM conventions, and any
   prior logs. Preserve user changes and existing outputs.
2. Required local inputs, transferred by the user into `Data/inputs/`, are:
   - `all_orthologs.csv`
   - `all_disease_orthologs.csv`
3. Use the tracked canonical gene ID tables:
   - `Data/outputs/sexspecific_disease_ids.csv`
   - `Data/outputs/all_disease_ortholog_ids.csv`
   - `Data/outputs/all_ortholog_ids.csv`
   - `Data/outputs/dmel_gene_reference.csv`
4. The local orthology CSVs contain FlyBase IDs, HGNC IDs, and DIOPT scores.
   The current data-wrangling code does not apply an explicit DIOPT threshold.
   Do not describe the set as filtered to a high-confidence DIOPT threshold
   unless the input provenance establishes that. Report the observed DIOPT
   score distributions and all pair-to-gene mapping attrition.
5. Public references to inspect and use if accessible:
   - Human protein FASTA, gene annotation, and xref from the same directory:
     `https://ftp.ebi.ac.uk/pub/ensemblorganisms/GCA/000/001/405/29/ensembl/2026_04/geneset/`
     The directory contains `pep.fa.bgz`, `genes.gtf.bgz` (or `.gz`), and
     `xref.tsv.gz`. `pep.fa.bgz` is BGZF-compressed; verify the exact file names
     and checksums from `md5sum.txt`.
   - Fly protein FASTA and annotation matching BDGP6.54 / Ensembl release 115:
     `https://ftp.ensembl.org/pub/release-115/fasta/drosophila_melanogaster/pep/Drosophila_melanogaster.BDGP6.54.pep.all.fa.gz`
     `https://ftp.ensembl.org/pub/release-115/gtf/drosophila_melanogaster/Drosophila_melanogaster.BDGP6.54.115.gtf.gz`
     Verify checksums and ensure FBgn and protein/transcript identifiers map
     unambiguously. These release files can replace a missing ignored local
     GTF; do not transfer the large local GTF unless the user asks.
6. Confirm the human xref/annotation actually maps the input HGNC IDs to genes
   and proteins and that fly annotation maps FBgn IDs to proteins. Record
   versions, retrieval dates, checksums, and unresolved/ambiguous mappings.
   Do not assume HPC network access; if downloads are blocked, report exact
   user-run transfer commands and stop before substituting another release.

## Analysis decisions and safeguards

- The existing canonical categories are mutually exclusive and prioritized:
  SSDOs; then disease IDs excluding SSDOs; then orthologue IDs excluding all
  disease IDs. Reproduce these sets from the tracked canonical ID files.
- Derive sequence comparisons only from the existing fly–human orthology
  pairs. Protein similarity must not be used to select or redefine orthologues.
- Reduce to one human orthologue per fly gene using the highest DIOPT score;
  resolve equal-score ties deterministically by sorted stable HGNC ID, without
  using identity, coverage, or the outcome values. Preserve the full pair
  table and document this rule.
- Use a consistent representative-protein rule for both species: Ensembl
  canonical protein where a reliable mapping is available; otherwise the
  longest valid protein for that gene, with stable protein ID as tie-break.
  Report missing and ambiguous mappings. Also provide a sensitivity analysis
  restricted to unambiguous one-to-one fly–human gene pairs.
- Prefer a documented pairwise global protein alignment method available in
  the HPC environment. BLASTP is a local aligner: if used, do not call its
  `pident` a whole-protein identity without qualification. Consolidate
  overlapping HSPs without double-counting residues; define identity and
  fly/human coverage explicitly and retain HSP/alignment details for audit.
  Do not silently install software or switch reference releases.
- Produce one observation per fly gene per metric. Report eligible genes,
  selected pairs, alignments, exclusions, and missingness by category. Identity
  and coverage have metric-specific measured universes; do not force them to
  share one universe or compare across different universes without disclosure.
- Plot group bootstrap means with 95% percentile intervals (5,000 resamples,
  seed 1) for identity and coverage. Include sample sizes and a clear axis/unit.
- For each metric separately, run the three two-sided pairwise differences in
  arithmetic means among the three categories using 10,000 two-sample
  permutations (seeded reproducibly). Apply BH correction across the three
  contrasts separately within identity and within coverage. Export raw and
  adjusted P-values, group sizes, mean differences, permutation count, seed,
  and method. Do not interpret non-significance as equivalence.
- Make intermediate pair-level and gene-level tables auditable. Keep bulk
  sequences, alignments, logs, and large intermediates untracked; document
  inputs, outputs, commands, definitions, and exclusions in a concise
  side-directory README.

## Execution protocol

- Do lightweight file/module/reference checks on the login node only.
- Use `ncbi-blast/2.16.0+` only if the selected alignment method requires it;
  inspect available modules and existing SLURM examples before choosing
  resources. Do not infer account, partition, paths, or resource limits.
- Before submitting any job, present the exact `sbatch` command, script,
  required inputs, expected outputs, scientific assumptions, and uncertainties
  to the user. Do not submit until the user explicitly approves that exact
  command.
- After approval, submit once, confirm the job ID and initial pending/running
  status, then stop monitoring. Do not inspect running logs, submit dependent
  stages, or overwrite results absent a later explicit request.
- Do not commit generated data or credentials. At completion, report exact
  output paths and provide explicit user-run `scp` commands for the compact
  gene-level tables, summary/test tables, plots, and README needed for local
  review. The user will pull/rebase before transferring outputs.
