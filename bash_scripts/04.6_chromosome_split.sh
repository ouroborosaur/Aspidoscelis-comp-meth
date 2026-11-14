#!/bin/bash
#SBATCH -J 04.6_chromosome_split
#SBATCH -p htc
#SBATCH -q public
#SBATCH -t 0-01:00:00
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
PLOT_DIR="${OUTDIR}/plot_ready"
LOGS_DIR="${OUTDIR}/logs"
mkdir -p "$PLOT_DIR" "$LOGS_DIR"
# input file
FINAL="${PLOT_DIR}/${SPECIES}_final.bed"

# confirm input file exists and is not empty
if [ ! -s "$FINAL" ]; then
    echo "ERROR: Final BED file not found or empty at ${FINAL}."
    exit 1
fi

# activate environment
module load mamba/latest
source activate $METHYLATION_ENV

echo "Splitting ${SPECIES} final BED by chromosome on $(date)..."

CHROM=$(cut -f1 "${PLOT_DIR}/${SPECIES}_final.bed" | tail -n +2 | sort | uniq)

# confirm at least one chromosome name was found
if [ -z "$CHROM" ]; then
    echo "ERROR: No chromosome names found in ${FINAL}."
    exit 1
fi

echo "Chromosomes detected: $(echo "$CHROM" | wc -l)"

# write each chromosome to separate BED file
for c in ${CHROM}; do
    grep "^${c}\b" "$FINAL" > "${PLOT_DIR}/${SPECIES}_${c}.bed"
    echo "  ${c}: $(wc -l < ${PLOT_DIR}/${SPECIES}_${c}.bed) sites"

done

source deactivate

echo "${SPECIES} chromosome splitting complete on $(date)."
