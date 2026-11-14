#!/bin/bash
#SBATCH -J 04.5_methylation_filter
#SBATCH -p htc
#SBATCH -q public
#SBATCH -t 0-00:30:00
#SBATCH -c 2
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
OUTDIR="${PROCESSED_DIR}/${SPECIES}"
METHYL_DIR="${OUTDIR}/methylation"
PLOT_DIR="${OUTDIR}/plot_ready"
LOGS_DIR="${OUTDIR}/logs"
mkdir -p "$METHYL_DIR" "$PLOT_DIR" "$LOGS_DIR"

# activate environment
module load mamba/latest
source activate $METHYLATION_ENV

# read in coverage cutoff file
CUTOFF_FILE="${OUTDIR}/logs/${SPECIES}_coverage_cutoff.txt"

# return error message if file doesn't exist
if [ ! -f "$CUTOFF_FILE" ]; then
    echo "ERROR: Coverage cutoff file not found at ${CUTOFF_FILE}."
    exit 1
fi

# read in coverage cutoff value
CUTOFF=$(cat "$CUTOFF_FILE")

# return error message if value is empty
if [ -z "$CUTOFF" ]; then
    echo "ERROR: Coverage cutoff value is empty."
    exit 1
fi

echo "Using coverage cutoff (mean + 3*stdev) ${CUTOFF} on $(date)..."

# input, temp files, final file
INPUT="${METHYL_DIR}/${SPECIES}_CpG.bed"
TMP1="${PLOT_DIR}/${SPECIES}_methylated.bed"
TMP2="${PLOT_DIR}/${SPECIES}_cov_filtered.bed"
FINAL="${PLOT_DIR}/${SPECIES}_final.bed"

# confirm input file exists and is not empty
if [ ! -s "$INPUT" ]; then
    echo "ERROR: Input CpG BED file not found or empty at ${INPUT}."
    exit 1
fi

echo "Filtering methylation calls for ${SPECIES} on $(date)..."

# keep only 5mC calls (modification code "m", excludes 5hmC "h")
awk '$4 == "m"' "$INPUT" > "$TMP1"
# print row count
echo "Rows after 5mC filter: $(wc -l < $TMP1)"

# filter for minimum coverage = 10, maximum coverage = dynamic cutoff
awk -v cutoff="$CUTOFF" '$5 >= 10 && $5 <= cutoff' "$TMP1" > "$TMP2"
# print row count
echo "Rows after coverage filter (>= 10 and <= ${CUTOFF}): $(wc -l < $TMP2)"

# reformat to tab-delimited, select columns
awk -v OFS='\t' '{print $1,$2,$3,$4,$6,$11,$5}' "$TMP2" > "$FINAL"

# add header
sed -i '1i chrom\tstart\tend\tmethyl\tstrand\tprob\tcov' "$FINAL"

# remove temp files
rm "$TMP1" "$TMP2"

# confirm final file is not empty
if [ ! -s "$FINAL" ]; then
    echo "ERROR: Final BED file is empty after filtering. Check cutoff value and input."
    exit 1
fi

echo "Final site count (excluding header): $(tail -n +2 $FINAL | wc -l)"

source deactivate

echo "${SPECIES} methylation filter complete on $(date)."
