#!/usr/bin/env bash
#SBATCH --job-name=dada2_16S
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=72:00:00
#SBATCH --mem=120G
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

set -euo pipefail

cd "$(git rev-parse --show-toplevel)"
mkdir -p logs

set +u
export MAMBA_ROOT_PREFIX="$PWD/.mamba"
source <("$PWD/bin/micromamba" shell hook -s bash)
micromamba activate "$PWD/envs/r16s"
set -u

# Force R to use your user library where you installed Bioc 3.20 + dada2 from source
export R_LIBS_USER="$HOME/R/bioc_r44"
export R_LIBS_SITE=""

which Rscript
Rscript --version

Rscript scripts/R/02_Build_Physeq_Objects.R
