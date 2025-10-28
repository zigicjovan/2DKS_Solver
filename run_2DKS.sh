#!/bin/bash
#SBATCH --account=def-bprotas                # RAC account for Bartek, to be used for priority jobs
#SBATCH --mail-user=zigicj@mcmaster.ca       # email address
#SBATCH --mail-type=ALL                      # Send email regarding all events with job
#SBATCH --ntasks=1                           # Serial code - keep at 1
#SBATCH --cpus-per-task=8                    
#SBATCH --mem=0M		                     # total memory required; default unit is megabytes
#SBATCH --time=03-00:00		                 # time (DD-HH:MM) or (HH:MM:SS)
#SBATCH --job-name=K0ell112120               # keep track of jobs

module load matlab/2024b.1

matlab -nodisplay -nodesktop -nosplash -r "main_2DKS(1e-03,24,0,0,1,1.12,1.20,0.02,1.0,2.0,11); exit" > K0ell112120.log 2>&1 < /dev/null
