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
K_start =           4.0
K_end =             4.0
K_step =            0.5
num_modes_start =   1.98#np.sqrt(13)#
num_modes_end =     num_modes_start#np.sqrt(20)#
ell1_start =        num_modes_start
ell1_end =          num_modes_end
ell1_step =         0.12
ell2_start =        ell1_start
ell2_end =          ell1_end
ell2_step =         ell1_step
K_range =           np.round(np.arange(K_start, K_end + K_step/2, K_step), 1)
ell1_range =        np.round(np.arange(ell1_start, ell1_end + ell1_step/2, ell1_step), 2)
#ell2_range =        np.round(np.arange(ell2_start, ell2_end + ell2_step/2, ell2_step), 2)
#K_range = np.round(np.array([3.5]), 1)
#ell1_range = [round(x, 2) for x in [np.sqrt(2),np.sqrt(3),np.sqrt(4),np.sqrt(6),np.sqrt(8),np.sqrt(9),np.sqrt(10),np.sqrt(13),np.sqrt(16), 
#               np.sqrt(17),np.sqrt(18),np.sqrt(19),np.sqrt(20),np.sqrt(23),np.sqrt(26),np.sqrt(29),np.sqrt(32)]]
#ell1_range = [round(x, 2) for x in [np.sqrt(13),np.sqrt(16),np.sqrt(18),np.sqrt(20)]]
ell2_range = ell1_range

# User-editable global settings
SBATCH_TIME = "03-00:00"  # requested time limit (D-HH:MM) 
RUN_ARRAY_NAME = f"b_{K_start}_{K_end}-{ell1_start:.2f}_{ell1_end:.2f}_{ell2_start:.2f}_{ell2_end:.2f}-branch.sh"

def generate_tasks():
    # Parameter choice formulas for 2DKS problem
    N_choice =  [ 48,   64,   96,   128,  160,  192,  256,  320,  384,  512,  576,  648,  720,  768,  810,  864,  900,  972,  1024 ]
    dt_choice = [ 5e-4, 2e-4, 1e-4, 5e-5, 2e-5, 1e-5, 5e-6, 2e-6, 1e-6, 5e-7, 4e-7, 3e-7, 2e-7, 1e-7, 9e-8, 8e-8, 7e-8, 6e-8, 5e-8 ]
    K_ref = 3.5
    T_width = round(0.4, 2)
    T_step = round(0.02, 2)

    # generate parameter tuples
    tasks = []
    filecount = 0.0
    memcount = 0.0
    for K in K_range:
        for ell1 in ell1_range:
            #for ell2 in ell2_range:
            ell2 = ell1
            elltemp = ell1
            T_target = 0.5*ell1 + 2.1
            #targettemp = 1.02*np.sqrt(1)
            #T_target = 0.5*targettemp + 2.1
            #T_range = np.round(np.array([(T_target / 2) - 3*K ]), 2) # symmetry
            #T_range = np.round(np.array([(2*T_target / 3) - K ]), 2) # init
            #T_range = np.round(np.arange(0.7, 1.4, 0.1), 2) # fixed
            #T_range = np.round(np.array([-2.85,-2.83,-2.81,-2.79,-2.77,-2.75,-2.71,-2.69,-2.67,-2.65,-2.63,-2.61,-2.59,-2.57,-2.55]), 2) # fixed
            #T_range = np.round(np.array([(T_target - T_width) - K ]), 2) # Lower Bound LB
            #T_range = np.round(np.array([(T_target + T_width) - K ]), 2) # Upper Bound UB
            T_range = np.round(np.arange((T_target - T_width) - K, (T_target + T_width) - K + T_step/2, T_step), 2) # [LB,UB] branch
            idx = max( 0 , min( int( np.round(elltemp + 2*(K - K_ref) + 3) ) , len(N_choice) - 1) ) 
            N = N_choice[idx] 
            dt = dt_choice[idx]
            for T in T_range:
                mem_est = 8 * np.ceil((1.1 * np.exp( -4.789714989 + 0.83721882 * np.log10((N**2) * (10.0**T) / dt) - 1.70503490 * elltemp
                                            + 0.23432402 * K + 0.18691756 * elltemp * np.log10((N**2) * (10.0**T) / dt) )) / 8.0 )
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
        array_spec = f"0-{count-1}%1"  # for sequential tasks
        #array_spec = f"0-{count-1}"  # for concurrent tasks

        lines.append(f"echo \"Group memory={mem}: tasks={count}, param_file={param_file}\"")
        lines.append("if [[ $DRY_RUN -eq 0 ]]; then")
        lines.append("  if [[ -z \"$prev_jobid\" ]]; then")
        lines.append(f"    prev_jobid=$(sbatch --parsable --account=def-bprotas --mail-user=zigicj@mcmaster.ca --mail-type=ALL \\")
        lines.append(f"                     --job-name={RUN_ARRAY_NAME} \\")
        lines.append(f"                     --ntasks=1 --cpus-per-task=8 --time={SBATCH_TIME} --mem={mem} \\")
        lines.append(f"                     --array={array_spec} --output=slurm_logs/slurm-%A_%a.out --error=slurm_logs/slurm-%A_%a.err \\")
        lines.append(f"                     --export=ALL,PARAM_FILE={param_file} ./run_task_array.sh)")
        lines.append("  else")
        lines.append(f"    prev_jobid=$(sbatch --parsable --dependency=afterany:${{prev_jobid}} --account=def-bprotas --mail-user=zigicj@mcmaster.ca --mail-type=ALL \\")
        lines.append(f"                     --job-name={RUN_ARRAY_NAME} \\")
        lines.append(f"                     --ntasks=1 --cpus-per-task=8 --time={SBATCH_TIME} --mem={mem} \\")
        lines.append(f"                     --array={array_spec} --output=slurm_logs/slurm-%A_%a.out --error=slurm_logs/slurm-%A_%a.err \\")
        lines.append(f"                     --export=ALL,PARAM_FILE={param_file} ./run_task_array.sh)")
        lines.append("  fi")
        lines.append("  echo \"Submitted job: $prev_jobid\"")
        lines.append("else")
        lines.append("  if [[ -z \"$prev_jobid\" ]]; then")
        lines.append(f"    echo \"Dry run: would sbatch --parsable --job-name={RUN_ARRAY_NAME} --mem={mem} --cpus-per-task=8 --time={SBATCH_TIME} --array={array_spec} --export=ALL,PARAM_FILE={param_file} ./run_task_array.sh\"")
        lines.append("  else")
        lines.append(f"    echo \"Dry run: would sbatch --parsable --dependency=afterany:$prev_jobid --job-name={RUN_ARRAY_NAME} --mem={mem} --cpus-per-task=8 --time={SBATCH_TIME} --array={array_spec} --export=ALL,PARAM_FILE={param_file} ./run_task_array.sh\"")
        lines.append("  fi")
        lines.append("fi")
        lines.append("")

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
