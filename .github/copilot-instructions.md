# Cross-project lab instructions

- Minimal explanations. Bullets over paragraphs. 
- Work conservatively on scientific code and analyses. 
- Do not change core analytical logic, statistical definitions, grouping schemes, or biological interpretation without explicit approval.
- Preserve provenance carefully. When data are merged across studies, keep input sources, matched universes, and filtering logic explicit.
- For genomics and transcriptomics workflows, prefer exporting intermediate tables that make downstream figures and tests auditable.
- When reporting statistical results, state the relevant universe, numerator, denominator, and any study-matching or exclusion logic when those affect interpretation.
- Prefer existing project conventions for Fisher's exact enrichment tests, two-sample permutation tests, bootstrap means, and confidence intervals unless those conventions are statistically inappropriate for the question or data structure; in that case, explain the problem clearly and propose or implement the least disruptive correct fix.
- If a result depends on a specific universe or assay background, do not compare groups across mismatched universes without explaining the consequence.
- Avoid introducing new languages, frameworks, or file formats into an established repo unless clearly justified by the existing project style.
- Code may be edited locally or directly in an HPC clone. Do not assume HPC access from a local session; when running on the HPC, use SLURM for substantial computation and follow the submission protocol below.
- Prioritize reproducibility, explicit paths/inputs, and clear annotations over clever automation.
- In research repos, keep README and inline comments accurate, minimal, and reviewer-friendly.
- If a request is ambiguous but a conservative interpretation is possible, choose the least disruptive scientifically sensible option and state that choice clearly.

## Local and HPC data handling

- Check Git state before editing or running workflows. A clone may lack ignored inputs, large outputs, exploratory files, and uncommitted changes present on another machine.
- Preserve the repository's existing tracked figures and summary-table conventions. Keep large reference downloads, alignment files, and other bulk HPC inputs/outputs untracked unless explicitly requested otherwise; add new bulk-data paths to `.gitignore` when introducing a workflow.
- Use explicit repository-relative input/output paths. Document the specific files to transfer between HPC and local clones when downstream analyses require them.
- Do not automate, capture, request, or expose passwords. For transfers, provide a user-run `scp` command with explicit source and destination paths; the user authenticates interactively.

## HPC execution and troubleshooting

- Keep workflows runnable and reviewable by scientists. Normal submission should be `sbatch path/to/script.sh` from the repository root.
- Put resource requests, job names, and log paths in `#SBATCH` headers. Put module loads, input/output paths, and scientific settings in a clearly labelled configuration section rather than long submission-time environment overrides. Ensure log directories exist before submission.
- Follow existing SLURM conventions. Prefer a short, annotated analysis-specific entry script calling shared code over duplicated logic or a new configuration framework. Document purpose, inputs, outputs, and the one-line command succinctly in the relevant README.
- Keep exploratory workflows separate from production. Do not promote them or change scientific defaults without approval.
- Before changing an HPC workflow or submitting a job, check required inputs, reference/annotation releases and identifier compatibility, expected outputs, modules/executables, and relevant existing SLURM stdout/stderr. Do not infer machine-specific paths, accounts, partitions, or resource limits from another project.
- Run substantial computations through SLURM, not on login nodes. Request a valid partition when necessary; do not select a compute node unless explicitly directed by the user.
- Do not delete, overwrite, or regenerate existing HPC results without explicit approval. Fail explicitly with clear diagnostics rather than guessing fixes.
- Retain failed-job logs during diagnosis. After a successful replacement and output verification, remove only superseded failed-job logs with user approval. Keep successful-job logs.
- Record durable execution findings, versions, commands, and scientifically relevant failure causes in the appropriate tracked documentation. Do not commit credentials or machine-private information.

## SLURM submission protocol

- Before submission, present the exact command, required inputs and paths, expected outputs, scientific assumptions, and remaining uncertainty. Submit only after the user explicitly approves that specific command.
- After an approved submission, confirm once that SLURM accepted the job and that it is pending or running, using `squeue` and, when useful, `sacct`. Report the job ID and initial state, then stop monitoring and finish the submission task. Job acceptance is not evidence that the analysis succeeded.
- The user will return to request inspection after completion or failure. Do not continuously poll jobs, inspect running-job logs, submit dependent stages, overwrite results, or regenerate outputs unless the user explicitly requests that follow-up work.