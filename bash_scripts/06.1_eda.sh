#!/bin/bash
#SBATCH -J 06.1_eda
#SBATCH -p htc
#SBATCH -q public
#SBATCH -t 0-04:00:00
#SBATCH -c 4
#SBATCH --mem=64G
#SBATCH -o /scratch/maelstr1/thesis/logs/%x_%j.out
#SBATCH -e /scratch/maelstr1/thesis/logs/%x_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=maelstr1@asu.edu
#SBATCH --export=NONE

# exit immediately on command failure (-e), unset variables (-u), or failure
#  anywhere in a pipeline
set -euo pipefail

source /scratch/maelstr1/thesis/scripts/config.sh

# optional first positional argument:
# FALSE (default) = use files containing all retained sites
# TRUE = use repeat-filtered files where they exist
# repeat-filtered angusticeps is currently unavailable
# TRUE run contains deppii & cozumelus data only
# cozumelus_ang is retained for descriptive EDA only
REPEAT_FILTERED="${1:-FALSE}"

# reject misspelled/unsupported values before starting R job
case "${REPEAT_FILTERED^^}" in
    TRUE|FALSE) ;;
    *)
        echo "ERROR: argument must be TRUE or FALSE"
        exit 1
        ;;
esac

# mkdir -p is safe if LOG_DIR already exists; it does not erase its contents.
mkdir -p "$LOG_DIR"

echo "Coordinate-independent EDA starting on $(date)..."
echo "Repeat filtered: ${REPEAT_FILTERED}"

module load mamba/latest
source activate "$METHYLKIT_ENV"


Rscript "${ANALYSIS_DIR}/r_scripts/06.1_eda.R" \
    --processed_dir "$PROCESSED_DIR" \
    --output_dir "${ANALYSIS_DIR}/results" \
    --repeat_filtered "$REPEAT_FILTERED" \
    2>&1 | tee "${LOG_DIR}/eda_${REPEAT_FILTERED}_$(date +%Y%m%d).log"


source deactivate

echo "Coordinate-independent EDA completed on $(date)."
