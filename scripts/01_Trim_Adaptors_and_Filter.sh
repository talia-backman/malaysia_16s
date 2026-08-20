#!/usr/bin/env bash
#SBATCH --job-name=cutadapt_16S
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=24:00:00
#SBATCH --mem=64G
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
export R_LIBS_USER="$HOME/R/bioc_r44"
export R_LIBS_SITE=""

which Rscript
Rscript --version

Rscript scripts/R/01_trim_filter_reads.R
