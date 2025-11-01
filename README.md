# Animations of output with various parameter settings:
https://www.youtube.com/playlist?list=PLwsovxEJkjzJrUHeRQMrPbvsQRfkz0e-W

# Main function to run branch of tests for various K, ell, T values:
main_2DKS(dtc,Nc,Kstart,Kend,Knum,ellstart,ellend,ellgap,Tstart,Tend,Tnum,runc,continuationc,tolc,optTc)

# To run quick test on Linux terminal with disabled respawning (assuming 'output' folder exists):
nohup matlab -nodisplay -nodesktop -nosplash -r "main_2DKS(.0001,32,0,0,1,1.02,1.02,0.02,0,0,1); exit" > ./output/output1.log 2>&1 < /dev/null &

# To run on HPC cluster:
sbatch run_2DKS.sh

# To get PID from queue and check status, or to cancel:
sq
seff PID
scancel PID

# To fetch new code update from Github (but leave all untracked files):
git fetch origin
git reset --hard origin/main

# Check storage:
du -h -d 3 /scratch/zigicj/2DKS_Solver/ | sort -h
diskusage_report

# Copy to scratch for production runs, and set up tests there:
cp -r /project/def-bprotas/zigicj/2DKS_Solver/* /scratch/zigicj/2DKS_Solver/
python3 param_driver.py
bash run_array.sh --dry-run
bash run_array.sh

# Monitor output file in terminal:
tail -f ./output/output1.log

# List all running process PIDs:
ps aux | grep matlab | grep -v grep

# Kill all active processes:
pkill -f "main_2DKS"

# If necessary, force kill respawning active process:
# Step 1 - check current PID parent process
pstree -ps PID
# Step 2 - identify:
systemd(...)
 └─ systemd(...)
     └─ MATLAB(parentPID)
# Step 3 - kill parent process:
kill -9 parentPID

# Push update to github:
./gitpush.sh . "" "COMMIT MESSAGE"

