#!/bin/bash
#SBATCH -J 04.1_quality_filter
#SBATCH -p htc
#SBATCH -q public
#SBATCH -t 0-1:30:00
#SBATCH -c 20
#SBATCH --mem=100G
#SBATCH -o /scratch/maelstr1/thesis/logs/%x_%j.out
#SBATCH -e /scratch/maelstr1/thesis/logs/%x_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user="maelstr1@asu.edu"
#SBATCH --export=NONE

source /scratch/maelstr1/thesis/scripts/config.sh

# pass species as argument
SPECIES=$1
ID=$2

# directories
BAM="${PROCESSED_DIR}/${SPECIES}_${ID}/aligned/${SPECIES}_${ID}_sorted.bam"
REF="${REF_DIR}/fasta/${SPECIES}_ref.fasta"
OUTDIR="${PROCESSED_DIR}/${SPECIES}_${ID}"
LOGS_DIR="${OUTDIR}/logs"
FILTERED_BAM="${OUTDIR}/aligned/${SPECIES}_${ID}_filtered.bam"
mkdir -p "$LOGS_DIR"

# activate environment
module load samtools-1.16-gcc-11.2.0
module load mamba/latest
source activate $METHYLATION_ENV


# QUALITY FILTER
echo "Filtering BAM for ${SPECIES}_${ID} on $(date)..."

samtools view \
    -@ $SLURM_CPUS_PER_TASK \
    -b -F 2308 -q 20 \
    "$BAM" \
    -o "$FILTERED_BAM"


# INDEX
echo "Indexing ${SPECIES} on $(date)..."

samtools index "$FILTERED_BAM"

echo "Indexing complete on $(date)."


source deactivate

echo "${SPECIES}_${ID} BAM filtering complete on $(date)."
