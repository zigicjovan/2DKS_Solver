### Animations of output with various parameter settings:
https://www.youtube.com/playlist?list=PLwsovxEJkjzJrUHeRQMrPbvsQRfkz0e-W

### To run quick test on Linux terminal (assuming 'output' folder exists):
1. make run

### To run on HPC cluster:
1. Modify param_driver.py
2. [Terminal] python3 param_driver.py
3. Run execution script from terminal output of param_driver.py

### To generate figures on cluster:
1. cd Data
2. module load matlab/2024b.1
3. matlab -nodisplay -nosplash
4. testDir = '[folder name]';
5. generateFigures(testDir);

### To get PID from queue to check status or to kill process:
1. sq
2. scancel PID

### Check storage and priority:
1. du -h -d 3 /user/foldername/ | sort -h
2. diskusage_report
3. sshare -l -A def-prof1_cpu -u prof1,grad2,postdoc3
