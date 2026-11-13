#!/bin/bash
set -euo pipefail
cd "$(dirname "$0")/.."
PYTHON="fly_human_conservation/.venv_local/bin/python"
[[ -x "$PYTHON" ]] || { echo "Create the local Python environment first: $PYTHON" >&2; exit 1; }
printf -v WORKFLOW_WRAPPER_COMMAND '%q ' "$0" "$@"
export WORKFLOW_WRAPPER_COMMAND PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
exec "$PYTHON" -u fly_human_conservation/align_proteins.py "$@"
