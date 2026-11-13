#!/bin/bash
set -euo pipefail
cd "$(dirname "$0")/.."
PYTHON="fly_human_conservation/.venv_local/bin/python"
[[ -x "$PYTHON" ]] || { echo "Create the local Python environment first: $PYTHON" >&2; exit 1; }
printf -v WORKFLOW_WRAPPER_COMMAND '%q ' "$0" "$@"
export WORKFLOW_WRAPPER_COMMAND PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1
exec "$PYTHON" -u fly_human_conservation/prepare_references.py \
    --apply-reviewed-recoveries fly_human_conservation/recovery_review_4955170/tables/per_gene_resolution.csv "$@"
