#!/bin/bash
#SBATCH --job-name=fly_human_recovered_align
#SBATCH --partition=compute
#SBATCH --account=bio
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=08:00:00
#SBATCH --output=fly_human_conservation/logs/recovered-align-%j.out
#SBATCH --error=fly_human_conservation/logs/recovered-align-%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=rpr23sxu@uea.ac.uk

# Submit from the repository root after approval:
# sbatch fly_human_conservation/align_recovered.sh
set -euo pipefail

module load python/anaconda/2024.10/3.12.7
module -t list 2>&1
export PYTHONNOUSERSITE=1
export PYTHONDONTWRITEBYTECODE=1
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1

PROJECT_DIR="${SLURM_SUBMIT_DIR:-$PWD}"
SIDE="fly_human_conservation"
MAPPING="$SIDE/mapping_4955446"
ENVIRONMENT="side/fly_human_conservation/.venv_alignment_biopython186_numpy1264"

cd "$PROJECT_DIR"
if [[ -z "${SLURM_JOB_ID:-}" ]]; then
    echo "Error: submit this alignment through sbatch, not on the login node." >&2
    exit 1
fi
if [[ ! -s "$MAPPING/provenance/PREPARATION_COMPLETE.json" ]]; then
    echo "Error: missing completed recovery mapping; submit from the repository root." >&2
    exit 1
fi
if [[ ! -x "$ENVIRONMENT/bin/python" ]]; then
    echo "Error: the pinned Biopython environment is missing; do not install automatically." >&2
    exit 1
fi
if [[ -e "$SIDE/alignment_${SLURM_JOB_ID}" ]]; then
    echo "Error: fresh alignment output already exists; refusing to overwrite." >&2
    exit 1
fi

"$ENVIRONMENT/bin/python" -m pip --isolated check
"$ENVIRONMENT/bin/python" -m unittest discover -s "$SIDE" -p test_align_proteins.py -v
"$ENVIRONMENT/bin/python" -u "$SIDE/align_proteins.py" --mapping-dir "$MAPPING"
