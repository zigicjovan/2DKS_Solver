#!/bin/bash
DRY_RUN=0
if [[ "$1" == "--dry-run" ]]; then DRY_RUN=1; fi
mkdir -p slurm_logs output runscripts
echo "(run_array.sh) dry-run=$DRY_RUN"

echo "Group memory=32G: tasks=84, param_file=runscripts/task_params_32G.txt"
if [[ $DRY_RUN -eq 0 ]]; then
  sbatch --account=def-bprotas --mail-user=zigicj@mcmaster.ca --mail-type=ALL \
         --ntasks=1 --cpus-per-task=8 --time=07-00:00 --mem=32G \
         --array=0-83 --output=slurm_logs/slurm-%A_%a.out --error=slurm_logs/slurm-%A_%a.err \
         --export=ALL,PARAM_FILE=runscripts/task_params_32G.txt ./run_task_array.sh
else
  echo "Dry run: would sbatch --array=0-83 --mem=32G --cpus-per-task=8 --export=PARAM_FILE=runscripts/task_params_32G.txt ./run_task_array.sh"
fi
