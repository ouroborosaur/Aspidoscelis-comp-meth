#!/bin/bash
#SBATCH -J build_methylation_env
#SBATCH -p htc
#SBATCH -q public
#SBATCH -t 0-02:00:00
#SBATCH -N 1
#SBATCH -c 2
#SBATCH --mem=8G
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
ENV_PATH="/scratch/maelstr1/conda_envs/methylation"

# remove old environment
echo "Removing old methylation environment on $(date)..."
mamba env remove -p "$ENV_PATH" -y 2>/dev/null || true
# clean up directory
rm -rf "$ENV_PATH"

# create environment with required packages
echo "Creating methylation environment on $(date)..."
mamba create -p "$ENV_PATH" \
    -c bioconda \
    -c conda-forge \
    samtools \
    pod5 \
    ont-modkit \
    bedtools \
    bcftools \
    -y

# confirm environment directory was created
if [ ! -d "${ENV_PATH}/bin" ]; then
    echo "ERROR: Environment directory not created at ${ENV_PATH}."
    echo "Installation failed."
    exit 1
fi

echo "Installation complete on $(date). Verifying..."

# verify environment building worked
source activate "$ENV_PATH"

# verify each package with explicit exit on failure to avoid silent passing
echo "Checking samtools..."
if ! samtools --version | head -1; then
    echo "ERROR: samtools not found after installation."
    exit 1
fi

echo "Checking modkit..."
if ! modkit --version; then
    echo "ERROR: modkit not found after installation."
    exit 1
fi

echo "Checking bedtools..."
if ! bedtools --version | head -1; then
    echo "ERROR: bedtools not found after installation."
    exit 1
fi

echo "Checking Python path..."
if ! python -c "import encodings; print('Python path valid')"; then
    echo "ERROR: Python standard library not accessible."
    exit 1
fi

source deactivate

echo "methylation environment setup complete on $(date)."
echo "Environment path: ${ENV_PATH}"
