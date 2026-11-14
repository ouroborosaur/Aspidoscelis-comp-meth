#!/bin/bash
#SBATCH -J 02_ref_prep
#SBATCH -p htc
#SBATCH -q public
#SBATCH -t 0-01:00:00
#SBATCH -c 4
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
REF_GZ="${REF_DIR}/fasta/${SPECIES}_ref.fasta.gz"
REF_FA="${REF_DIR}/fasta/${SPECIES}_ref.fasta"

# activate environment
module load mamba/latest
source activate $METHYLATION_ENV

echo "Decompressing reference for ${SPECIES} at $(date)..."
gunzip -k "$REF_GZ"

echo "Indexing with samtools faidx at $(date)..."
samtools faidx "$REF_FA"

source deactivate

echo "Reference ${REF_FA} ready at $(date)."
