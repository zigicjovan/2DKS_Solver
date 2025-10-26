# Animations of output with various parameter settings:
https://www.youtube.com/playlist?list=PLwsovxEJkjzJrUHeRQMrPbvsQRfkz0e-W

# Main function to run branch of tests for various K, ell, T values (check main_2DKS function for argument functionality):
main_2DKS(dtc,Nc,Kstart,Kend,Knum,ellstart,ellend,ellgap,Tstart,Tend,Tnum)

# Disable sleep (when running on local machine):
sudo systemctl mask sleep.target suspend.target hibernate.target hybrid-sleep.target

# Re-enable sleep (when running on local machine):
sudo systemctl unmask sleep.target suspend.target hibernate.target

# To run quick test on Linux terminal with disabled respawning:
export MW_NO_SERVICE_HOST=1
export MW_DISABLE_SERVICE_HOST=1
nohup matlab -nodisplay -nodesktop -nosplash -r "main_2DKS(.0001,32,0,0,1,1.02,1.02,0.02,0,0,1); exit" > output1.log 2>&1 < /dev/null &

# Monitor output file in terminal:
tail -f output1.log

# List all running process PIDs:
ps aux | grep matlab | grep -v grep

# If respawning disabled, force kill active process:
kill -9 PID

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
