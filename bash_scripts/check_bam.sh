#!/bin/bash
#SBATCH -J bam_readcount
#SBATCH -p htc
#SBATCH -q public
#SBATCH -t 0-00:45:00
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -o /scratch/maelstr1/thesis/logs/%x_%j.out
#SBATCH -e /scratch/maelstr1/thesis/logs/%x_%j.err
#SBATCH --export=NONE

source /scratch/maelstr1/thesis/scripts/config.sh

# pass species as argument
SPECIES=$1

# files
SORTED="${PROCESSED_DIR}/${SPECIES}/aligned/${SPECIES}_sorted.bam"
FILTERED="${PROCESSED_DIR}/${SPECIES}/aligned/${SPECIES}_filtered.bam"

# activate environment
module load mamba/latest
source activate methylation

echo "Sorted BAM read count:"
samtools view -c -@ 4 "$SORTED"

echo "Filtered BAM read count:"
samtools view -c -@ 4 "$FILTERED"

echo "Sorted BAM flagstat:"
samtools flagstat -@  4 "$SORTED"

echo "Filtered BAM flagstat:"
samtools flagstat -@ 4 "$FILTERED"

source deactivate
