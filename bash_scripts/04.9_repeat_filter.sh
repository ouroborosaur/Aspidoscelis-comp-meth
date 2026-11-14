#!/bin/bash
#SBATCH -J 04.9_repeat_filter
#SBATCH -p htc
#SBATCH -q public
#SBATCH -t 0-02:00:00
#SBATCH -c 4
#SBATCH --mem=32G
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
OUTDIR="${PROCESSED_DIR}/${SPECIES}_${ID}"
PLOT_DIR="${OUTDIR}/plot_ready"
LOGS_DIR="${OUTDIR}/logs"
FINAL="${PLOT_DIR}/${SPECIES}_${ID}_final.bed"
FINAL_NOREPEAT="${PLOT_DIR}/${SPECIES}_${ID}_final_norepeat.bed"
REPEATS_BED="${REPEATS_DIR}/${SPECIES}_repeats.bed"

mkdir -p "$PLOT_DIR" "$LOGS_DIR"

# confirm input BED exists
if [ ! -s "$FINAL" ]; then
    echo "ERROR: final.bed not found."
    exit 1
fi

# confirm repeats BED exists
if [ ! -f "$REPEATS_BED" ]; then
    echo "ERROR: RepeatMasker BED file not found."
    exit 1
fi

# activate environment
module load mamba/latest
source activate $METHYLATION_ENV

echo "Filtering repetitive regions for ${SPECIES}_${ID} on $(date)..."

# copy header
head -1 "$FINAL" > "$FINAL_NOREPEAT"

# bedtools intersect -v excludes sites overlapping RepeatMasker region
# -a input sites
# -b repeat regions to exclude
# -v keep sites in -a that do not overlap -b
tail -n +2 "$FINAL" | \
    bedtools intersect \
    -v \
    -a stdin \
    -b "$REPEATS_BED" \
    >> "$FINAL_NOREPEAT"

# report how many sites were removed
BEFORE=$(tail -n +2 "$FINAL" | wc -l)
AFTER=$(tail -n +2 "$FINAL_NOREPEAT" | wc -l)
REMOVED=$((BEFORE - AFTER))
PCT=$(echo "scale=1; $REMOVED * 100 / $BEFORE" | bc)

echo "  Sites before repeat filtering: ${BEFORE}"
echo "  Sites after repeat filtering: ${AFTER}"
echo "  Sites removed: ${REMOVED} (${PCT}%)"

# confirm output is not empty
if [ ! -s "$FINAL_NOREPEAT" ]; then
    echo "ERROR: final_norepeat.bed is empty after filtering."
    exit 1
fi

# subgenome-specific for cozumelus
if [ "$SPECIES" == "cozumelus" ]; then
    echo "Creating subgenome-specific norepeat files for cozumelus on $(date)..."
    ANG_SUBG_DIR="${PLOT_DIR}/subgenome_angusticeps"
    DEP_SUBG_DIR="${PLOT_DIR}/subgenome_deppii"
    mkdir -p "$ANG_SUBG_DIR" "$DEP_SUBG_DIR"

    # angusticeps norepeat file
    head -1 "$FINAL_NOREPEAT" > "${ANG_SUBG_DIR}/cozumelus_ang_subgenome_norepeat.bed"
    grep "^angusticeps_Chr_" "$FINAL_NOREPEAT" >> "${ANG_SUBG_DIR}/cozumelus_ang_subgenome_norepeat.bed"

    ANG_NOREPEAT=$(tail -n +2 "${ANG_SUBG_DIR}/cozumelus_ang_subgenome_norepeat.bed" | wc -l)
    echo "  A. angusticeps subgenome norepeat sites: ${ANG_NOREPEAT}"

    # deppii norepeat file
    head -1 "$FINAL_NOREPEAT" > "${DEP_SUBG_DIR}/cozumelus_dep_subgenome_norepeat.bed"
    grep "^deppii_Chr_" "$FINAL_NOREPEAT" >> "${DEP_SUBG_DIR}/cozumelus_dep_subgenome_norepeat.bed"

    DEP_NOREPEAT=$(tail -n +2 "${DEP_SUBG_DIR}/cozumelus_dep_subgenome_norepeat.bed" | wc -l)
    echo "  A. deppii subgenome norepeat sites: ${DEP_NOREPEAT}"
fi

source deactivate

echo "Repeat filtering for ${SPECIES}_${ID} complete on $(date)."
