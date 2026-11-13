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
- Prioritize reproducibility, explicit paths/inputs, and clear annotations over clever automation.
- In research repos, keep README and inline comments accurate, minimal, and reviewer-friendly.
- If a request is ambiguous but a conservative interpretation is possible, choose the least disruptive scientifically sensible option and state that choice clearly.

## Data handling

- Check Git state before editing or running workflows. A clone may lack ignored inputs, large outputs, exploratory files, and uncommitted changes present on another machine.
- Preserve the repository's existing tracked figures and summary-table conventions. Keep large reference downloads, alignment files, and other bulk inputs/outputs untracked unless explicitly requested otherwise; add new bulk-data paths to `.gitignore` when introducing a workflow.
- Use explicit repository-relative input/output paths. Keep workflows runnable and reviewable by scientists, and document their purpose, inputs, outputs and commands succinctly.
- Keep exploratory workflows separate from production. Do not promote them or change scientific defaults without approval.
- Before running workflows, check required inputs, reference/annotation releases and identifier compatibility, expected outputs, and required executables.
- Do not delete, overwrite, or regenerate existing results without explicit approval. Fail explicitly with clear diagnostics rather than guessing fixes.
- Record durable execution findings, versions, commands, and scientifically relevant failure causes in the appropriate tracked documentation. Do not commit credentials or machine-private information.