#!/bin/bash
#SBATCH --account=def-bprotas                # RAC account for Bartek, to be used for priority jobs
#SBATCH --mail-user=zigicj@mcmaster.ca       # email address
#SBATCH --mail-type=ALL                      # Send email regarding all events with job
#SBATCH --ntasks=1                           # Serial code - keep at 1
#SBATCH --cpus-per-task=8                    
#SBATCH --mem=32G		             # total memory required; default unit is megabytes
#SBATCH --time=05-00:00		             # time (DD-HH:MM) or (HH:MM:SS)
#SBATCH --job-name=K5ell106-26               # keep track of jobs

module load matlab/2024b.1

matlab -nodisplay -nodesktop -nosplash -r "main_2DKS(5e-08,192,5,5,1,1.06,1.06,0.02,-2.6,-2.6,1,'optimize','IC',1e-6,0.0); exit" > ./output/K5ell106-26.log 2>&1 < /dev/null
