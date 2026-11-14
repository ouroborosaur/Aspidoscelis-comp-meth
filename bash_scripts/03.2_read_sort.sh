#!/bin/bash
#SBATCH -J 03.2_read_sort
#SBATCH -p public
#SBATCH -q public
#SBATCH -t 0-06:00:00
#SBATCH -c 16
#SBATCH --mem=64G
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
ALIGNED_BAM="${OUTDIR}/aligned/${SPECIES}_aligned.bam"
SORTED_BAM="${OUTDIR}/aligned/${SPECIES}_sorted.bam"
mkdir -p "${OUTDIR}/aligned" "${OUTDIR}/logs"

# activate environment
module load htslib-1.21-gcc-12.1.0
module load mamba/latest
source activate $METHYLATION_ENV


# SORT
echo "Sorting ${SPECIES} on $(date)..."

samtools sort \
    -@ $SLURM_CPUS_PER_TASK \
    -o "$SORTED_BAM" \
    "$ALIGNED_BAM"

echo "Sorting completed on $(date)."


# INDEX
echo "Indexing ${SPECIES} on $(date)..."

samtools index "$SORTED_BAM"

echo "Indexing complete on $(date)."


# QUALITY CHECK
# flagstat counts reads by alignment flag category
# check primary mapped percentage (target >85%)
# check secondary & supplementary counts (non-zero before filtering)
samtools flagstat \
    "$SORTED_BAM" \
    > "${OUTDIR}/logs/flagstat.txt"

source deactivate

echo "Sorting & indexing for ${SPECIES} complete on $(date)."
