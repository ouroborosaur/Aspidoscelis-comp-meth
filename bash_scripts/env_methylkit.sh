#!/bin/bash
#SBATCH -J build_methylkit_env
#SBATCH -p htc
#SBATCH -q public
#SBATCH -t 0-04:00:00
#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -o /scratch/maelstr1/thesis/logs/%x_%j.out
#SBATCH -e /scratch/maelstr1/thesis/logs/%x_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user="maelstr1@asu.edu"
#SBATCH --export=NONE

# clear path conflicts
unset PYTHONPATH
unset PYTHONHOME

export CONDA_PKGS_DIRS="/scratch/maelstr1/conda_pkgs"
mkdir -p /scratch/maelstr1/conda_pkgs

module load mamba/latest

# verify install to SCRATCH, not HOME
PKGS_DIR=$(conda config --show pkgs_dirs 2>/dev/null | grep scratch | head -1)
if [ -z "$PKGS_DIR" ]; then
    echo "ERROR: conda pkgs_dirs is not configured to use scratch."
    exit 1
fi

# define installation path
ENV_PATH="/scratch/maelstr1/conda_envs/methylkit_env"

# remove old environment
echo "Removing old methylkit environment on $(date)..."
mamba env remove -p "$ENV_PATH" -y 2>/dev/null || true
# clean up directory
rm -rf "$ENV_PATH"

# create environment with required packages
echo "Creating methylkit_env environment on $(date)..."
mamba create -p "$ENV_PATH" \
    -c bioconda \
    -c conda-forge \
    r-base \
    bioconductor-methylkit \
    bioconductor-genomicranges \
    r-optparse \
    -y

# confirm environment directory was created
if [ ! -d "${ENV_PATH}/bin" ]; then
    echo "ERROR: Environment directory not created at ${ENV_PATH}."
    echo "Installation failed."
    exit 1
fi

echo "Installation complete on $(date). Verifying..."

source activate "$ENV_PATH"

echo "Checking R libraries..."
if ! Rscript -e "
    library(methylKit)
    library(GenomicRanges)
    library(optparse)
    cat('methylKit version:', as.character(packageVersion('methylKit')), '\n')
    cat('R environment valid\n')
"; then
    echo "ERROR: R library check failed."
    exit 1
fi

source deactivate

echo "methylkit_env setup complete on $(date)."
echo "Environment path: ${ENV_PATH}"
