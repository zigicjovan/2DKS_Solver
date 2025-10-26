# Animations of output with various parameter settings:
# https://www.youtube.com/playlist?list=PLwsovxEJkjzJrUHeRQMrPbvsQRfkz0e-W

main_2DKS(dt,N,range of K values,range of ell values,range of T values)

# Disable sleep
sudo systemctl mask sleep.target suspend.target hibernate.target hybrid-sleep.target

# Re-enable sleep
sudo systemctl unmask sleep.target suspend.target hibernate.target

# To run quick test on linux terminal:
nohup matlab -nodisplay -nodesktop -nosplash -r "main_2DKS(.0001,32,0,0,1,1.02,1.02,0.02,0,0,1); exit" > output1.log 2>&1 < /dev/null &

# Monitor output file in terminal:
tail -f output1.log

# List all running processes:
ps aux | grep matlab

# If necessary, kill or force kill the process:
pkill -f "main_2DKS"
pkill -9 -f "main_2DKS"

# Push update to github:
./gitpush.sh . "" "COMMIT MESSAGE"
