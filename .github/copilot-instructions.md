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
- For HPC-backed projects, assume code editing happens locally, SLURM/HPC execution done manually by user, large HPC outputs may be untracked and copied manually to the local respective directories. Do not assume direct access to protected HPC environments.
- Prioritize reproducibility, explicit paths/inputs, and clear annotations over clever automation.
- In research repos, keep README and inline comments accurate, minimal, and reviewer-friendly.
- If a request is ambiguous but a conservative interpretation is possible, choose the least disruptive scientifically sensible option and state that choice clearly.