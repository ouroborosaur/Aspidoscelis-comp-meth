#!/bin/bash
#SBATCH -J 06.2_stats
#SBATCH -p htc
#SBATCH -q public
#SBATCH -t 0-02:00:00
#SBATCH -c 4
#SBATCH --mem=32G
#SBATCH -o /scratch/maelstr1/thesis/logs/%x_%j.out
#SBATCH -e /scratch/maelstr1/thesis/logs/%x_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=maelstr1@asu.edu
#SBATCH --export=NONE

set -euo pipefail
source /scratch/maelstr1/thesis/scripts/config.sh

# optional first argument
# FALSE (default) = all retained sites
# TRUE = repeat filtered

REPEAT_FILTERED="${1:-FALSE}"
case "${REPEAT_FILTERED^^}" in
    TRUE|FALSE) ;;
    *)
        echo "ERROR: argument must be TRUE or FALSE"
        exit 1
        ;;
esac

mkdir -p "$LOG_DIR"

echo "Coordinate-independent statistics starting on $(date)..."
echo "Repeat filtered: ${REPEAT_FILTERED}"

module load mamba/latest
source activate "$METHYLKIT_ENV"


Rscript "${ANALYSIS_DIR}/r_scripts/06.2_stats.R" \
    --output_dir "${ANALYSIS_DIR}/results" \
    --repeat_filtered "$REPEAT_FILTERED" \
    --bootstrap_reps 5000 \
    --pattern_threshold 10 \
    2>&1 | tee "${LOG_DIR}/stats_${REPEAT_FILTERED}_$(date +%Y%m%d).log"

source deactivate

echo "Coordinate-independent statistics completed on $(date)."
