#!/bin/bash
#SBATCH -J 05.1_methylartist_segmeth
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
ID=$2
# pass chromosome as optional argument (test with a single chromosome before running on full genome)
TEST_CHROM=${3:-""}
WINDOW=100000
# optional confidence thresholds (empty = methylartist defaults)
# approximate modkit's confidence thresholds (methylated, canonical)
THRESH_ARGS="--meth_thresh 0.75 --can_thresh 0.75"

OUTDIR="${PROCESSED_DIR}/${SPECIES}_${ID}"
FILTERED_BAM="${OUTDIR}/aligned/${SPECIES}_${ID}_filtered.bam"
REF="${REF_DIR}/fasta/${SPECIES}_ref.fasta"
MA_DIR="${OUTDIR}/methylartist"
mkdir -p "$MA_DIR"

# test runs do not overwrite full runs
if [ -n "$TEST_CHROM" ]; then
    SEGMETH_OUT="${MA_DIR}/${SPECIES}_${ID}_${TEST_CHROM}_test.segmeth.tsv"
else
    SEGMETH_OUT="${MA_DIR}/${SPECIES}_${ID}.segmeth.tsv"
fi

# exit on missing inputs
for f in "$FILTERED_BAM" "${FILTERED_BAM}.bai" "${REF}.fai"; do
    if [ ! -s "$f" ]; then
        echo "ERROR: Required file missing or empty: $f"
        exit 1
    fi
done

# build windows from reference index
# column 1 = name; column 2 = length; column 4 = label (subgenome for cozumelus, "genome" for sexual species)
SEG_BED="${MA_DIR}/${SPECIES}_${ID}_windows_${WINDOW}.bed"
awk -v w="$WINDOW" -v OFS='\t' '{
    if ($1 ~ /^angusticeps_Chr_/) lab = "ang_subgenome";
    else if ($1 ~ /^deppii_Chr_/) lab = "dep_subgenome";
    else if ($1 ~ /^Chr_/) lab = "genome";
    else lab = "unplaced";
    for (s = 0; s < $2; s += w) {
        e = (s + w < $2) ? s + w : $2;
        print $1, s, e, lab, ".";
    }
}' "${REF}.fai" > "$SEG_BED"

# optional single-chromosome test
if [ -n "$TEST_CHROM" ]; then
    awk -v c="$TEST_CHROM" '$1 == c' "$SEG_BED" > "$SEG_BED.tmp" && mv "${SEG_BED}.tmp" "$SEG_BED"
    echo "TESTING: restricted to ${TEST_CHROM}"
fi

if [ ! -s "$SEG_BED" ]; then
    echo "ERROR: window BED is empty. Check chromosome name against ${REF}.fai"
    exit 1
fi
echo "Windows to process: $(wc -l < $SEG_BED)"

unset PYTHONPATH
unset PYTHONHOME
module load mamba/latest
source activate $METHYLARTIST_ENV

# run from output directory
cd "$MA_DIR"
echo "segmeth for ${SPECIES}_${ID} starting on $(date)..."
methylartist segmeth \
    -b "$FILTERED_BAM" \
    -i "$SEG_BED" \
    -o "$SEGMETH_OUT" \
    -p $SLURM_CPUS_PER_TASK \
    -m m \
    --ref "$REF" \
    --motif CG \
    $THRESH_ARGS

# confirm a segmeth table was produced and has more than a header
if [ ! -s "$SEGMETH_OUT" ] || [ "$(wc -l < "$SEGMETH_OUT")" -lt 2 ]; then
    echo "ERROR: segmeth output missing or empty: ${SEGMETH_OUT}."
    exit 1
fi

source deactivate

echo "segmeth for ${SPECIES}_${ID} completed on $(date). Output: ${SEGMETH_OUT}"
