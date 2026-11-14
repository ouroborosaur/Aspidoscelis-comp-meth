#!/bin/bash
#SBATCH -J 05.2_methylartist_plots
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
ID=$2
# pass chromosome as optional (test) argument
TEST_CHROM=${3:-""}
MINCALLS=100

MA_DIR="${PROCESSED_DIR}/${SPECIES}_${ID}/methylartist"

# test runs do not overwrite full runs
if [ -n "$TEST_CHROM" ]; then
    SEGMETH_TSV="${MA_DIR}/${SPECIES}_${ID}_${TEST_CHROM}_test.segmeth.tsv"
else
    SEGMETH_TSV="${MA_DIR}/${SPECIES}_${ID}.segmeth.tsv"
fi

# exit on missing file
if [ ! -s "$SEGMETH_TSV" ]; then
    echo "ERROR: ${SEGMETH_TSV} not found."
    exit 1
fi

# for hybrid, only plot two subgenomes (omit unplaced contigs)
CAT_ARGS=""
if  [ "$SPECIES" == "cozumelus" ]; then
    CAT_ARGS="-c ang_subgenome,dep_subgenome"
else
    CAT_ARGS="-c genome"
fi

unset PYTHONPATH
unset PYTHONHOME
export MPLBACKEND=Agg    # draw straight to file
module load mamba/latest
source activate $METHYLARTIST_ENV

# run from output directory
cd "$MA_DIR"

echo "Plotting from $SEGMETH_TSV on $(date)..."

# strip plot
methylartist segplot -s "$SEGMETH_TSV" -n "$MINCALLS" $CAT_ARGS
# violin plot
methylartist segplot -s "$SEGMETH_TSV" -n "$MINCALLS" $CAT_ARGS -v

source deactivate

echo "Plots written on $(date) to ${MA_DIR}:"
ls -lh "$MA_DIR"/*.png
