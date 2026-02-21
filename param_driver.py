#!/usr/bin/env python3
"""
param_driver.py
Generates runscripts/task_params_<mem>_<timestamp>.txt and placeholder scripts,
and writes run_array.sh which submits arrays per-memory-group.
"""

import numpy as np
import os
from pathlib import Path
import time

# User-editable parameter ranges
K_start =  4.0
K_end =  4.0
K_step = 0.5
ell1_start = 1.02*np.sqrt(256)
ell1_end = ell1_start
ell1_step = 0.12
ell2_start = 1.02*np.sqrt(4)
ell2_end = ell2_start
ell2_step = 0.12
K_range = np.round(np.arange(K_start, K_end + K_step/2, K_step), 1)
ell1_range = np.round(np.arange(ell1_start, ell1_end + ell1_step/2, ell1_step), 2)
ell2_range = np.round(np.arange(ell2_start, ell2_end + ell2_step/2, ell2_step), 2)
#K_range = np.round(np.array([3.5]), 1)
#ell1_range = [round(x, 2) for x in [2.04]] 

# User-editable global settings
SBATCH_TIME = "02-00:00"  # requested time limit (D-HH:MM) 
RUN_ARRAY_NAME = f"{K_start}_{K_end}-{ell1_start:.2f}_{ell1_end:.2f}_{ell2_start:.2f}_{ell2_end:.2f}-branch.sh"

def generate_tasks():
    # Parameter choice formulas for 2DKS problem
    N_choice =  [ 48,   64,   96,   128,  160,  192,  256,  320,  384,  512,  576,  648,  720,  768,  810,  864,  900,  972,  1024 ]
    dt_choice = [ 5e-4, 2e-4, 1e-4, 5e-5, 2e-5, 1e-5, 5e-6, 2e-6, 1e-6, 5e-7, 4e-7, 3e-7, 2e-7, 1e-7, 9e-8, 8e-8, 7e-8, 6e-8, 5e-8 ]
    K_ref = 3.5
    T_width = round(0.02, 2)
    T_step = round(0.02, 2)

    # generate parameter tuples
    tasks = []
    filecount = 0.0
    memcount = 0.0
    for K in K_range:
        for ell1 in ell1_range:
            for ell2 in ell2_range:
                #T_target = 0.5*ell1 + 2.08
                targettemp = 1.8*np.sqrt(1)
                T_target = 0.5*targettemp + 2.08
                #T_range = np.round(np.array([(T_target / 2) - 3*K ]), 2) # init
                #T_range = np.round(np.array([(T_target - T_width) - K ]), 2) # Lower Bound LB
                T_range = np.round(np.array([(T_target + T_width) - K ]), 2) # Upper Bound UB
                #T_range = np.round(np.arange((T_target - T_width) - K, (T_target + T_width) - K + T_step/2, T_step), 2) # [LB,UB] branch
                elltemp = 1.02
                idx = max( 0 , min( int( np.round(elltemp + 2*(K - K_ref) + 3) ) , len(N_choice) - 1) ) 
                N = N_choice[idx] 
                dt = dt_choice[idx]
                for T in T_range:
                    mem_est = 8 * np.ceil((1.25 * np.exp( -4.789714989 + 0.83721882 * np.log10((N**2) * (10.0**T) / dt) - 1.70503490 * elltemp
                                                + 0.23432402 * K + 0.18691756 * elltemp * np.log10((N**2) * (10.0**T) / dt) ) + 2.0) / 8.0 )
                    mem_req = int(np.maximum( 32 , mem_est ))
                    mem = f"{mem_req}G" 
                    addmemcount = mem_req
                    addfilecount = 0.1 * 1e-3 * np.power(10.0,T) / dt
                    memcount += addmemcount
                    filecount += addfilecount
                    tasks.append((float(K), float(ell1), float(ell2), float(T), float(dt), int(N), mem))
    return tasks, int(filecount), int(memcount)

# Write runscripts/ and param files
def write_runscripts(tasks, out_dir='runscripts'):
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)

    # timestamp tag to avoid overwriting old param files
    tag = time.strftime("%Y%m%d_%H%M%S")

    # group tasks by memory
    groups = {}
    for idx, t in enumerate(tasks):
        mem = t[-1]
        groups.setdefault(mem, []).append((idx, t))

    # write task_params_<mem>_<tag>.txt
    for mem, items in groups.items():
        mem_tag = mem.replace('G', 'G')
        param_fname = out / f"task_params_{mem_tag}_{tag}.txt"
        with param_fname.open('w') as pf:
            for (idx, (K, ell1, ell2, T, dt, N, mem)) in items:
                output_file = f"output/run_{K}_{ell1:.2f}_{ell2:.2f}_{T:.2f}_{dt:.0e}_{N}.mat"
                pf.write(f"{idx} {K} {ell1:.2f} {ell2:.2f} {T:.2f} {dt:.0e} {N} {mem} {output_file}\n")

    # placeholder scripts
    for idx, (K, ell1, ell2, T, dt, N, mem) in enumerate(tasks):
        fname = out / f"run_{K}_{ell1:.2f}_{ell2:.2f}_{T:.2f}_{dt}_{N}.sh"
        with fname.open('w') as fh:
            fh.write(f"""#!/bin/bash
#SBATCH --account=def-bprotas
#SBATCH --mail-user=zigicj@mcmaster.ca
#SBATCH --mail-type=ALL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem={mem}
#SBATCH --time={SBATCH_TIME}
module load matlab/2024b.1
# Placeholder; real execution uses ../run_task_array.sh
""")
    return groups, tag

# Write run_array.sh driver
def write_run_array_sh(groups, tag, out_fname='run_array.sh', max_concurrent=0):
    lines = [
        "#!/bin/bash",
        'DRY_RUN=0',
        'if [[ "$1" == "--dry-run" ]]; then DRY_RUN=1; fi',
        'mkdir -p slurm_logs output runscripts',
        'echo "(run_array.sh) dry-run=$DRY_RUN"',
        ''
    ]

    for mem, items in groups.items():
        mem_tag = mem.replace('G', 'G')
        param_file = f"runscripts/task_params_{mem_tag}_{tag}.txt"
        count = len(items)
        array_spec = f"0-{count-1}%{max_concurrent}" if max_concurrent > 0 else f"0-{count-1}"

        lines.append(f"echo \"Group memory={mem}: tasks={count}, param_file={param_file}\"")
        lines.append("if [[ $DRY_RUN -eq 0 ]]; then")
        lines.append(f"  sbatch --account=def-bprotas --mail-user=zigicj@mcmaster.ca --mail-type=ALL \\")
        lines.append(f"         --job-name={RUN_ARRAY_NAME} \\")
        lines.append(f"         --ntasks=1 --cpus-per-task=8 --time={SBATCH_TIME} --mem={mem} \\")
        lines.append(f"         --array={array_spec} --output=slurm_logs/slurm-%A_%a.out --error=slurm_logs/slurm-%A_%a.err \\") # for concurrent tasks
        #lines.append(f"         --array={array_spec}%1 --output=slurm_logs/slurm-%A_%a.out --error=slurm_logs/slurm-%A_%a.err \\") # for sequential tasks
        lines.append(f"         --export=ALL,PARAM_FILE={param_file} ./run_task_array.sh")
        lines.append("else")
        lines.append(f"  echo \"Dry run: would sbatch --array={array_spec}%1 --mem={mem} --cpus-per-task=8 --export=PARAM_FILE={param_file} ./run_task_array.sh\"")
        lines.append("fi\n")

    Path(out_fname).write_text("\n".join(lines))
    os.chmod(out_fname, 0o755)

# Main
def main():
    [tasks,filecount,memcount] = generate_tasks()
    print(f"Parameter combinations: {len(tasks)}")
    print(f"Filespace required: {(filecount)}K")
    print(f"Memory requested: {(memcount)}GB")
    print("File count check: find . -type f -printf '%h\n' | sort | uniq -c | sort -nr")
    print("Disk usage check: du -h | sort -hr")
    groups, tag = write_runscripts(tasks)
    write_run_array_sh(groups, tag, out_fname=RUN_ARRAY_NAME, max_concurrent=0)
    print(f"Wrote runscripts/ (task param files with tag {tag})")
    print(f"Execution script: bash {RUN_ARRAY_NAME}")
    print("Priority check: sshare -l -A def-bprotas_cpu -u fabianbl,zigicj,noahb")

if __name__ == "__main__":
    main()
