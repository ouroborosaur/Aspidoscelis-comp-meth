#!/bin/bash
#SBATCH -J 03.1_read_align
#SBATCH -p public
#SBATCH -q public
#SBATCH -t 2-00:00:00
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH --gres=gpu:a100:1
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
POD5_DIR="${RAW_DIR}/pod5/${SPECIES}_${ID}"
REF="${REF_DIR}/fasta/${SPECIES}_ref.fasta"
OUTDIR="${PROCESSED_DIR}/${SPECIES}_${ID}"
mkdir -p "${OUTDIR}/aligned" "${OUTDIR}/logs"

# activate environment
module load htslib-1.21-gcc-12.1.0
module load cuda-12.6.1-gcc-12.1.0
module load mamba/latest
source activate $METHYLATION_ENV

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK

# basecalling, methylation calling, alignment...
echo "Basecalling & aligning ${SPECIES}_${ID} on $(date)..."

# --nv and --bind allow dorado to see the GPU and read & write files
#  super accurate basecalling for 5-Methylcytosine and 5-Hydroxymethylcytosine methylation
apptainer exec --nv \
    --bind /scratch/maelstr1:/scratch/maelstr1 \
    ${DORADO_SIF} \
    dorado basecaller \
    --device cuda:0 \
    --recursive \
    sup,5mCG_5hmCG \
    "$POD5_DIR" \
    --reference "$REF" \
    2> "${OUTDIR}/logs/basecall_align.log" \
    > "${OUTDIR}/aligned/${SPECIES}_${ID}_aligned.bam"

source deactivate

echo "Read alignment for ${SPECIES} completed on $(date)."
