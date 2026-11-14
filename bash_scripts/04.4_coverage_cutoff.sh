#!/bin/bash
#SBATCH -J 04.4_coverage_cutoff
#SBATCH -p htc
#SBATCH -q public
#SBATCH -t 0-00:10:00
#SBATCH -c 2
#SBATCH --mem=16G
#SBATCH -o /scratch/maelstr1/thesis/logs/%x_%j.out
#SBATCH -e /scratch/maelstr1/thesis/logs/%x_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user="maelstr1@asu.edu"
#SBATCH --export=NONE

source /scratch/maelstr1/thesis/scripts/config.sh

# pass species as argument
SPECIES=$1

# directories
OUTDIR="${PROCESSED_DIR}/${SPECIES}"
METHYL_DIR="${OUTDIR}/methylation"
LOGS_DIR="${OUTDIR}/logs"
mkdir -p "$METHYL_DIR" "$LOGS_DIR"

# activate environment
module load mamba/latest
source activate $METHYLATION_ENV

# calculate dynamic coverage cutoff
echo "Calculating coverage cutoff for ${SPECIES} on $(date)..."

CUTOFF=$(awk 'NR>1 {sum+=$10; sumsq+=$10^2; n++} END {mean=sum/n; stdev=sqrt(sumsq/n - mean^2); cutoff=mean + 3*stdev; print cutoff}' "${METHYL_DIR}/${SPECIES}_CpG.bed")

# confirm cutoff variable/value is not empty
if [ -z "$CUTOFF" ]; then
    echo "ERROR: Coverage cutoff calculation failed. Check that ${METHYL_DIR}/${SPECIES}_CpG.bed exists."
    exit 1
fi

echo "Coverage cutoff (mean + 3*stdev) for ${SPECIES}: ${CUTOFF}"

# save cutoff value for use downstream
echo "${CUTOFF}" > "${LOGS_DIR}/${SPECIES}_coverage_cutoff.txt"

source deactivate

echo "${SPECIES} coverage cutoff complete on $(date)."
