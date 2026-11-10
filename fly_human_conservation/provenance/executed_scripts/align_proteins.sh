#!/bin/bash
#SBATCH --job-name=fly_human_align
#SBATCH --partition=compute
#SBATCH --account=bio
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=08:00:00
#SBATCH --output=side/fly_human_conservation/logs/align-%j.out
#SBATCH --error=side/fly_human_conservation/logs/align-%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=rpr23sxu@uea.ac.uk

# Submit from the repository root:
# sbatch side/fly_human_conservation/align_proteins.sh
set -euo pipefail

# Configuration: fixed successful preparation, releases, and scoring are in Python.
module load python/anaconda/2024.10/3.12.7
module -t list 2>&1
PROJECT_DIR="${SLURM_SUBMIT_DIR:-$PWD}"
SIDE="side/fly_human_conservation"
ENVIRONMENT="$SIDE/.venv_alignment_biopython186_numpy1264"
export PYTHONNOUSERSITE=1
export PYTHONDONTWRITEBYTECODE=1
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1

cd "$PROJECT_DIR"
if [[ -z "${SLURM_JOB_ID:-}" ]]; then
    echo "Error: submit this alignment through sbatch, not on the login node." >&2
    exit 1
fi
if [[ ! -s "$SIDE/mapping_4935338/provenance/PREPARATION_COMPLETE.json" ]]; then
    echo "Error: missing successful preparation; submit from the repository root." >&2
    exit 1
fi
if [[ ! -e "$ENVIRONMENT" ]]; then
    python -m venv "$ENVIRONMENT"
    "$ENVIRONMENT/bin/python" -m pip --isolated --disable-pip-version-check install \
        --only-binary=:all: --requirement "$SIDE/alignment_requirements.txt" \
        --report "$ENVIRONMENT/install_report.json"
fi
"$ENVIRONMENT/bin/python" -m pip --isolated check
"$ENVIRONMENT/bin/python" -m unittest discover \
    -s "$SIDE" -p test_align_proteins.py -v
"$ENVIRONMENT/bin/python" -u "$SIDE/align_proteins.py"
