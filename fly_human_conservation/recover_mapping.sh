#!/bin/bash
#SBATCH --job-name=fly_human_recover_map
#SBATCH --partition=compute
#SBATCH --account=bio
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=fly_human_conservation/logs/recover-map-%j.out
#SBATCH --error=fly_human_conservation/logs/recover-map-%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=rpr23sxu@uea.ac.uk

# Submit from the repository root after approval:
# sbatch fly_human_conservation/recover_mapping.sh
set -euo pipefail

module load python/anaconda/2024.10/3.12.7
module -t list 2>&1
export PYTHONNOUSERSITE=1
export PYTHONDONTWRITEBYTECODE=1

PROJECT_DIR="${SLURM_SUBMIT_DIR:-$PWD}"
ENTRY_POINT="fly_human_conservation/prepare_references.py"
REVIEW="fly_human_conservation/recovery_review_4955170/tables/per_gene_resolution.csv"

cd "$PROJECT_DIR"
if [[ -z "${SLURM_JOB_ID:-}" ]]; then
    echo "Error: submit this recovery mapping through sbatch, not on the login node." >&2
    exit 1
fi
if [[ ! -f "$ENTRY_POINT" || ! -f "$REVIEW" ]]; then
    echo "Error: submit from the repository root with the approved review present." >&2
    exit 1
fi
if [[ -e "fly_human_conservation/mapping_${SLURM_JOB_ID}" ]]; then
    echo "Error: fresh mapping output already exists; refusing to overwrite." >&2
    exit 1
fi

python fly_human_conservation/check_handoff.py
python -m unittest discover -s fly_human_conservation -p test_prepare_references.py -v
python -m unittest discover -s fly_human_conservation -p test_review_mapping_recovery.py -v
python -u "$ENTRY_POINT" \
    --reuse-verified-references \
    --apply-reviewed-recoveries "$REVIEW"
