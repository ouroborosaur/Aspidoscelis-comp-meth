#!/bin/bash
#SBATCH -J index_filtered_bam
#SBATCH -p htc
#SBATCH -q public
#SBATCH -t 0-01:00:00
#SBATCH -c 16
#SBATCH --mem=32G
#SBATCH -o /scratch/maelstr1/thesis/logs/%x_%j.out
#SBATCH -e /scratch/maelstr1/thesis/logs/%x_%j.err
#SBATCH --export=NONE

source /scratch/maelstr1/thesis/scripts/config.sh

# pass species as argument
SPECIES=$1

# activate environment
module load mamba/latest
source activate methylation

echo "Indexing filtered BAM on $(date)..."

samtools index -@ 16 \
    "/scratch/maelstr1/thesis/02_processed/${SPECIES}/aligned/${SPECIES}_filtered.bam"

echo "Indexing complete on $(date)."

source deactivate
