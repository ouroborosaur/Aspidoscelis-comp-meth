#!/bin/bash
#SBATCH -J 04.7_subgenome_split
#SBATCH -p htc
#SBATCH -q public
#SBATCH -t 0-01:00:00
#SBATCH -c 2
#SBATCH --mem=8G
#SBATCH -o /scratch/maelstr1/thesis/logs/%x_%j.out
#SBATCH -e /scratch/maelstr1/thesis/logs/%x_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user="maelstr1@asu.edu"
#SBATCH --export=NONE

source /scratch/maelstr1/thesis/scripts/config.sh

# pass species as argument
SPECIES=$1
ID=$2

# only run script for hybrid (cozumelus)
if [ "$SPECIES" != "cozumelus" ]; then
    echo "INFO: Subgenome splitting is only applicable to hybrid species."
    exit 0
fi

# directories
OUTDIR="${PROCESSED_DIR}/${SPECIES}_${ID}"
PLOT_DIR="${OUTDIR}/plot_ready"
LOGS_DIR="${OUTDIR}/logs"
FINAL="${PLOT_DIR}/${SPECIES}_${ID}_final.bed"

# output directories for subgenomes
ANG_SUBG_DIR="${PLOT_DIR}/subgenome_angusticeps"
DEP_SUBG_DIR="${PLOT_DIR}/subgenome_deppii"

mkdir -p "$ANG_SUBG_DIR" "$DEP_SUBG_DIR" "$LOGS_DIR"

# confirm input file exists and is not empty
if [ ! -s "$FINAL" ]; then
    echo "ERROR: Final BED file not found or empty at ${FINAL}."
    exit 1
fi

echo "Splitting cozumelus final.bed by subgenome on $(date)..."
echo "Using chromosome name prefix to identify subgenome origin."

# ANGUSTICEPS SUBGENOME
echo "angusticeps subgenome: chromosomes prefixed with 'angusticeps_Chr_'"

# copy header row from final.bed to the output file
head -1 "$FINAL" \
    > "${ANG_SUBG_DIR}/cozumelus_ang_subgenome.bed"

# grep all rows with angusticeps prefix
grep "^angusticeps_Chr_" "$FINAL" \
    >> "${ANG_SUBG_DIR}/cozumelus_ang_subgenome.bed"

# count sites
ANG_SITES=$(tail -n +2 "${ANG_SUBG_DIR}/cozumelus_ang_subgenome.bed" | wc -l)
echo "  angusticeps subgenome sites: ${ANG_SITES}"

# DEPPII SUBGENOME
echo "deppii subgenome: chromosomes prefixed with 'deppii_Chr_'"

# copy header row
head -1 "$FINAL" \
    > "${DEP_SUBG_DIR}/cozumelus_dep_subgenome.bed"

# grep all rows with deppii prefix
grep "^deppii_Chr_" "$FINAL" \
    >> "${DEP_SUBG_DIR}/cozumelus_dep_subgenome.bed"

# count sites
DEP_SITES=$(tail -n +2 "${DEP_SUBG_DIR}/cozumelus_dep_subgenome.bed" | wc -l)
echo "  deppii subgenome sites: ${DEP_SITES}"

# VALIDATION

# confirm subgenome files are not empty
if [ "$ANG_SITES" -eq 0 ]; then
    echo "ERROR: No sites found for angusticeps subgenome."
    exit 1
fi

if [ "$DEP_SITES" -eq 0 ]; then
    echo "ERROR: No sites found for deppii subgenome."
    exit 1
fi

# check all sites were assigned to a subgenome
TOTAL_SITES=$(tail -n +2 "$FINAL" | wc -l)
ASSIGNED=$((ANG_SITES + DEP_SITES))

echo "Total sites in final.bed: ${TOTAL_SITES}"
echo "Sites assigned to subgenomes: ${ASSIGNED}"

if [ "$TOTAL_SITES" -ne "$ASSIGNED" ]; then
    UNASSIGNED=$((TOTAL_SITES - ASSIGNED))
    echo "WARNING: ${UNASSIGNED} sites are excluded from both subgenome files."
    echo "These sites have an unrecognized chromosome prefix."
fi

echo "cozumelus subgenome splitting complete on $(date)."
