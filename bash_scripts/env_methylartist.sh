#!/bin/bash
#SBATCH -J build_methlyartist_env
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

module load mamba/latest

# verify install to SCRATCH, not HOME
if ! grep -q "/scratch/maelstr1/conda_pkgs" "$HOME/.condarc" 2>/dev/null; then
    echo "ERROR: ~/.condarc does not redirect pkgs_dirs to SCRATCH directory."
    exit 1
fi

# define installation path
ENV_PATH="/scratch/maelstr1/conda_envs/methylartist_env"

# remove old environment
echo "Removing old methylartist environment on $(date)..."
mamba env remove -p "$ENV_PATH" -y 2>/dev/null || true
rm -rf "${ENV_PATH:?}"

# create environment
echo "Creating methylartist environment on $(date)..."
mamba create -p "$ENV_PATH" \
    -c bioconda \
    -c conda-forge \
    methylartist \
    -y

# confirm environment directory was created
if [ ! -d "${ENV_PATH}/bin" ]; then
    echo "ERROR: Environment directory was not created at ${ENV_PATH}."
    echo "Installation failed."
    exit 1
fi

echo "Installation complete on $(date). Verifying..."

# verify environment building worked
source activate "$ENV_PATH"

# verify each package with explicit exit on failure
echo "Checking Python standard library..."
if ! python -c "import encodings; print('Python path OK')"; then
    echo "ERROR: Python standard library not accessible."
    exit 1
fi

echo "Checking key Python libraries..."
if ! python -c "import pysam, matplotlib, pandas, numpy; print('Library imports OK')"; then
    echo "ERROR: a required Python library failed to import."
    exit 1
fi

echo "Checking methylartist commands..."
for cmd in segmeth segplot locus region; do
    if ! methylartist "$cmd" --help > /dev/null 2>&1; then
        echo "ERROR: 'methylartist ${cmd} --help' failed."
        exit 1
    fi
    echo "  methylartist ${cmd}: OK"
done

# save segmeth and segplot option lists to confirm flags
HELP_DIR="/scratch/maelstr1/thesis/logs"
methylartist segmeth --help > "${HELP_DIR}/help_methylartist_segmeth.txt" 2>&1
methylartist segplot --help > "${HELP_DIR}/help_methylartist_segplot.txt" 2>&1
echo "Option lists saved to ${HELP_DIR}/help_methylartist_*.txt"

echo "Installed methylartist version:"
mamba list -p "$ENV_PATH" methylartist

source deactivate

echo "methylartist environment setup complete on $(date)."
echo "Environment path: ${ENV_PATH}"
