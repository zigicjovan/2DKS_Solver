#!/bin/bash
#SBATCH --account=mth260100
#SBATCH --partition=shared
#SBATCH --job-name=postprocess
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00-01:59
#SBATCH --output=slurm_logs/postprocess_%j.out
#SBATCH --error=output/postprocess_%j.err

module load python/3.9.5
unset PYTHONPATH
source "$PROJECT/x-jzigic/2DKS_Solver/.venv-postprocess/bin/activate"
export MPLBACKEND=Agg

cd "$SCRATCH/Data"
bash ./collectFigures.sh
