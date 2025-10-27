#!/bin/bash
#SBATCH --account=def-bprotas                # RAC account for Bartek, to be used for priority jobs
#SBATCH --mail-user=zigicj@mcmaster.ca       # email address
#SBATCH --mail-type=ALL                      # Send email regarding all events with job
#SBATCH --nodes=1		             # number of nodes (32 cpus per node)
#SBATCH --ntasks-per-node=32	             # number of MPI processes
#SBATCH --mem=0M		             # total memory required; default unit is megabytes
#SBATCH --time=01-00:00		             # time (DD-HH:MM) or (HH:MM:SS)
#SBATCH --job-name=K3ell112120               # keep track of jobs

module load matlab/2024b.1

# Run your MATLAB script (without GUI)
nohup matlab -nodisplay -nodesktop -nosplash -r "main_2DKS(5e-05,48,3,3,1,1.12,1.20,0.02,-0.5,0.2,8); exit" > outK3ell112120.log 2>&1 < /dev/null &
seff $SLURM_JOBID		# print short summary of cpu and memory efficiency (for more details: sacct -j <jobid>)




