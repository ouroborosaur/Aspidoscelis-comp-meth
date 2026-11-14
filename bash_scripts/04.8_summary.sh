#!/bin/bash
#SBATCH -J 04.8_summary
#SBATCH -p htc
#SBATCH -q public
#SBATCH -t 0-01:00:00
#SBATCH -c 20
#SBATCH --mem=40G
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
FILTERED_BAM="${OUTDIR}/aligned/${SPECIES}_filtered.bam"

mkdir -p "$METHYL_DIR" "$LOGS_DIR"

# confirm filtered BAM exists
if [ ! -s "$FILTERED_BAM" ]; then
    echo "ERROR: Filtered BAM not found."
    exit 1
fi

# confirm BAI index exists
if [ ! -f "${FILTERED_BAM}.bai" ]; then
    echo "ERROR: BAM index not found."
    exit 1
fi

# activate environment
module load mamba/latest
source activate $METHYLATION_ENV

echo "Running modkit summary for ${SPECIES} on $(date)..."

$MODKIT_BIN summary \
    --threads $SLURM_CPUS_PER_TASK \
    --log-filepath "${LOGS_DIR}/${SPECIES}_modkit_summary.log" \
    "$FILTERED_BAM" \
    > "${METHYL_DIR}/${SPECIES}_summary.txt"

# confirm summary was written
if [ ! -s "${METHYL_DIR}/${SPECIES}_summary.txt" ]; then
    echo "ERROR: Modkit summary output is empty."
    exit 1
fi

# print summary to stdout
echo "Summary output for ${SPECIES}: "
cat "${METHYL_DIR}/${SPECIES}_summary.txt"

source deactivate

echo "Modkit summary for ${SPECIES} complete on $(date)."
