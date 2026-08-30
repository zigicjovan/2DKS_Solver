#!/usr/bin/env bash
#SBATCH --account=def-bprotas
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --output=output/postprocess_%j.out
#SBATCH --error=output/postprocess_%j.err

module load StdEnv/2023 matlab

bash ./Data/collectFigures.sh
