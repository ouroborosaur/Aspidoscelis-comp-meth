#!/bin/bash
#SBATCH -J 01.2_extract_pod5_tar
#SBATCH -p public
#SBATCH -q public
#SBATCH -t 1-00:00:00
#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -o /scratch/maelstr1/thesis/logs/%x_%j.out
#SBATCH -e /scratch/maelstr1/thesis/logs/%x_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user="maelstr1@asu.edu"
#SBATCH --export=NONE

source /scratch/maelstr1/thesis/scripts/config.sh

# pass species and files as arguments
SPECIES=$1
TAR_DIR=$2

POD5_DEST="${RAW_DIR}/pod5/${SPECIES}"
TMP_EXTRACT="${RAW_DIR}/pod5/.tmp_extract_${SPECIES}"
mkdir -p "$POD5_DEST" "$TMP_EXTRACT"

# input validation
if [ -z "$SPECIES" ] || [ -z "$TAR_DIR" ]; then
    echo "ERROR: Missing required arguments."
    exit 1
fi

if [ ! -d "$TAR_DIR" ]; then
    echo "ERROR: Tar directory not found at ${TAR_DIR}."
    exit 1
fi

# find all tar's in the directory
mapfile -t TAR_FILES < <(find "$TAR_DIR" \
    -maxdepth 1 \
    -type f \
    \( -name "*.tar" -o -name "*.tar.gz" -o -name "*.tgz" \) \
    | sort)

# confirm at least one file was found
if [ ${#TAR_FILES[@]} -eq 0 ]; then
    echo "ERROR: No TAR files found in ${TAR_DIR}."
    echo "Contents of directory:"
    ls -lh "$TAR_DIR"
    exit 1
fi

# report what will be processed
echo "Species: ${SPECIES}"
echo "TAR directory: ${TAR_DIR}"
echo "TAR files found: ${#TAR_FILES[@]}"
echo "POD5 destination: ${POD5_DEST}"
echo ""

# extract tar's sequentially
TOTAL_POD5=0
FAILED_TARS=()
TAR_NUM=0

for TAR_FILE in "${TAR_FILES[@]}"; do
    TAR_NUM=$((TAR_NUM + 1))
    TAR_NAME=$(basename "$TAR_FILE")
    
    echo "  Processing TAR ${TAR_NUM} of ${#TAR_FILES[@]}: ${TAR_NAME}"
    echo "  Started: $(date)"
    
    # clean staging directory before each extraction
    rm -rf "${TMP_EXTRACT:?}"/*
    
    # determine extraction flags based on extension
    if [[ "$TAR_FILE" == *.tar.gz ]] || [[ "$TAR_FILE" == *.tgz ]]; then
        TAR_FLAGS="-xzf"
        echo "  Format: gzip-compressed TAR"
    else
        TAR_FLAGS="-xf"
        echo "  Format: uncompressed TAR"
    fi

    # extract to staging directory
    tar $TAR_FLAGS "$TAR_FILE" \
        -C "$TMP_EXTRACT" \
        --checkpoint=10000 \
        --checkpoint-action=echo="Extracted %u files from ${TAR_NAME}..."
    
    # check tar exit code; record tars for which extraction failed
    if [ $? -ne 0 ]; then
        echo "  ERROR: Extraction failed for ${TAR_NAME}."
        FAILED_TARS+=("$TAR_NAME")
        rm -rf "${TMP_EXTRACT:?}"/*
        continue
    fi

    # count and move pod5's to destination
    echo "Moving POD5 files to ${POD5_DEST} on $(date)..."
    POD5_THIS_TAR=0
    while IFS= read -r -d '' pod5_file; do
        DEST_FILE="${POD5_DEST}/$(basename "$pod5_file")"
        if [ -f "$DEST_FILE" ]; then
            TAR_PREFIX="${TAR_NAME%.tar.gz}"
            TAR_PREFIX="${TAR_PREFIX%.tgz}"
            TAR_PREFIX="${TAR_PREFIX%.tar}"
            DEST_FILE="${POD5_DEST}/${TAR_PREFIX}_$(basename "$pod5_file")"
            echo "  WARNING: Filename collision. Renamed to $(basename "$DEST_FILE")"
        fi
    mv "$pod5_file" "$DEST_FILE"
    POD5_THIS_TAR=$((POD5_THIS_TAR + 1))

    done < <(find "$TMP_EXTRACT" -type f -name "*.pod5" -print0)

    echo "  POD5 files extracted from ${TAR_NAME}: ${POD5_THIS_TAR}"
    echo "  Completed $(date)."
    echo ""

    TOTAL_POD5=$((TOTAL_POD5 + POD5_THIS_TAR))
    rm -rf "${TMP_EXTRACT:?}"/*

done

# final cleanup
rm -rf "$TMP_EXTRACT"

# summary report
echo "TAR EXTRACTION SUMMARY FOR ${SPECIES}"
echo "Tar files processed: ${TAR_NUM}"
echo "Total POD5 files extracted: ${TOTAL_POD5}"

if [ ${#FAILED_TARS[@]} -gt 0 ]; then
    echo ""
    echo "WARNING: The following TAR files failed to extract:"
    for FAILED in "${FAILED_TARS[@]}"; do
        echo "  - ${FAILED}"
    done
    echo "These archives may be corrupted."
    exit 1
fi

# confirm at least one pod5 landed in the destination
if [ "$TOTAL_POD5" -eq 0 ]; then
    echo "ERROR: No POD5 files were extracted."
    exit 1
fi

# report destination contents
echo ""
echo "POD5 files now in ${POD5_DEST}:"
ls "${POD5_DEST}"/*.pod5 2>/dev/null | wc -l | xargs -I{} echo "  {} files"
DEST_SIZE=$(du -sh "$POD5_DEST" | cut -f1)
echo "  Total size: ${DEST_SIZE}"

echo ""
echo "POD5 extraction for ${SPECIES} complete on $(date)."
