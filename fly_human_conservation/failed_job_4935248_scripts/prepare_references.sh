#!/bin/bash
#SBATCH --job-name=fly_human_prepare
#SBATCH --partition=compute
#SBATCH --account=bio
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=side/fly_human_conservation/logs/prepare-%j.out
#SBATCH --error=side/fly_human_conservation/logs/prepare-%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=rpr23sxu@uea.ac.uk

# Submit from the repository root:
# sbatch side/fly_human_conservation/prepare_references.sh
set -euo pipefail

# Configuration: fixed reference releases and paths are in the Python entry point.
module load python/anaconda/2024.10/3.12.7
module -t list 2>&1
export PYTHONNOUSERSITE=1
PROJECT_DIR="${SLURM_SUBMIT_DIR:-$PWD}"
ENTRY_POINT="side/fly_human_conservation/prepare_references.py"

cd "$PROJECT_DIR"
if [[ -z "${SLURM_JOB_ID:-}" ]]; then
    echo "Error: submit this preparation through sbatch, not on the login node." >&2
    exit 1
fi
if [[ ! -f "$ENTRY_POINT" ]]; then
    echo "Error: submit from the SSDOs_Howatson_Grieshop repository root." >&2
    exit 1
fi
python -u "$ENTRY_POINT"
