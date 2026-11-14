#!/bin/bash

# source /scratch/maelstr1/thesis/scripts/config.sh

# prevent path conflicts when activating conda environments
unset PYTHONPATH
unset PYTHONHOME

# environments
export CONDA_PKGS_DIRS="/scratch/maelstr1/conda_pkgs"
export METHYLATION_ENV="/scratch/maelstr1/conda_envs/methylation"
export METHYLKIT_ENV="/scratch/maelstr1/conda_envs/methylkit_env"
export METHYLARTIST_ENV="/scratch/maelstr1/conda_envs/methylartist_env"

# directories
export PROJECT_DIR="/scratch/maelstr1/thesis"
export SCRIPTS_DIR="${PROJECT_DIR}/scripts"
export RAW_DIR="${PROJECT_DIR}/00_raw"
export REF_DIR="${PROJECT_DIR}/01_references"
export REPEATS_DIR="${REF_DIR}/repeats"
export PROCESSED_DIR="${PROJECT_DIR}/02_processed"
export ANALYSIS_DIR="${PROJECT_DIR}/03_analysis"
export LOG_DIR="${PROJECT_DIR}/logs"

# tools
export DORADO_SIF="${PROJECT_DIR}/dorado_0.9.6.sif"
export DORADO_BIN="apptainer exec --nv --bind /scratch/maelstr1:/scratch/maelstr1 ${DORADO_SIF} dorado"
export MODKIT_BIN="/scratch/maelstr1/conda_envs/methylation/bin/modkit"
export DORADO_MODEL="${REF_DIR}/models/dna_r10.4.1_e8.2_400bps_sup@v5.0.0"
export CLAIR3_MODEL="${REF_DIR}/models/r1041_e82_400bps_sup_v500"

export THREADS=16
export MIN_COVERAGE=10
