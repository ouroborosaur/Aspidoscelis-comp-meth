#!/bin/bash
#SBATCH -J 04.3_pileup_cpg
#SBATCH -p htc
#SBATCH -q public
#SBATCH -t 0-00:30:00
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
REF="${REF_DIR}/fasta/${SPECIES}_ref.fasta"
OUTDIR="${PROCESSED_DIR}/${SPECIES}_${ID}"
METHYL_DIR="${OUTDIR}/methylation"
LOGS_DIR="${OUTDIR}/logs"
FILTERED_BAM="${OUTDIR}/aligned/${SPECIES}_${ID}_filtered.bam"
mkdir -p "$METHYL_DIR" "$LOGS_DIR"

# activate environment
module load samtools-1.16-gcc-11.2.0
module load mamba/latest
source activate $METHYLATION_ENV

echo "Running modkit pileup (CpG) for ${SPECIES}_${ID} on $(date)..."

$MODKIT_BIN pileup \
    "$FILTERED_BAM" \
    "${METHYL_DIR}/${SPECIES}_${ID}_CpG.bed" \
    --ref "$REF" \
    --modified-bases 5mC 5hmC \
    --mod-thresholds m:0.75 \
    --mod-thresholds h:0.75 \
    --filter-threshold 0.75 \
    --cpg \
    --threads $SLURM_CPUS_PER_TASK \
    2> "${LOGS_DIR}/modkit_cpg.log"

echo "Modkit pileup (CpG) complete on $(date)."

source deactivate
