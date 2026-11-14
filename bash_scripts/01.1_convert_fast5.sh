#!/bin/bash
#SBATCH -J 01.1_convert_fast5
#SBATCH -p public
#SBATCH -q public
#SBATCH -t 0-08:00:00
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH -o /scratch/maelstr1/thesis/logs/%x_%j.out
#SBATCH -e /scratch/maelstr1/thesis/logs/%x_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user="maelstr1@asu.edu"
#SBATCH --export=NONE

source /scratch/maelstr1/thesis/scripts/config.sh

# pass species as argument
SPECIES=$1

# directories
FAST5_IN="${RAW_DIR}/fast5/${SPECIES}"
POD5_OUT="${RAW_DIR}/pod5/${SPECIES}"
mkdir -p "$POD5_OUT"

# activate environment
module load mamba/latest
source activate $METHYLATION_ENV

echo "Converting FAST5 to POD5 for ${SPECIES} at $(date)... "
pod5 convert fast5 \
    "${FAST5_IN}"/*.fast5 \
    --output "${POD5_OUT}/reads.pod5" \
    --threads $SLURM_CPUS_PER_TASK

source deactivate

echo "Conversion completed at $(date)."

