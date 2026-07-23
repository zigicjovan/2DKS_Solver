#!/usr/bin/env python3
"""
param_driver.py
Generates runscripts/task_params_<mem>_<ranks>r_<timestamp>.txt and run scripts,
and writes an array driver grouped by memory and MPI rank count.
"""

import math
import numpy as np # type: ignore
import os
from pathlib import Path
import time

# User-editable parameter ranges
K_start =           4.0
K_end =             4.0
K_step =            0.5
num_modes_start =   1
num_modes_end =     num_modes_start
ell1_start =        1.35
ell1_end =          1.35
ell1_step =         0.15
ell2_start =        ell1_start
ell2_end =          ell1_end 
ell2_step =         0.1
K_range =           np.round(np.arange(K_start, K_end + K_step/2, K_step), 1)
ell1_range =        np.round(np.arange(ell1_start, ell1_end + ell1_step/2, ell1_step), 2)
# ell1_range = [round(x, 2) for x in [np.sqrt(13),np.sqrt(16),np.sqrt(18),np.sqrt(20)]]
ell2_range =        np.round(np.arange(ell2_start, ell2_end + ell2_step/2, ell2_step), 2)
# ell2_range = ell1_range

# User-editable global settings
SBATCH_TIME = "03-00:00"  # requested time limit (D-HH:MM) 
RUN_ARRAY_NAME = f"b_{K_start}_{K_end}-{ell1_start:.2f}_{ell1_end:.2f}_{ell2_start:.2f}_{ell2_end:.2f}-1.sh"
SEQUENTIAL_TASKS = False

# Solver settings shared by all generated production runs
IC = "s1"
OPTIMIZE = 1
TOL = "1e-10"
CONTINUATION = 0
OPT_T = 1
SAVE_STATES = 100

def generate_tasks():
    # Parameter choice formulas for 2DKS problem
    N_choice =  [ 64,   128,  192,  256,  320,  384,  512,  576,  640,  768,  896,  1024, 1536, 2048, 3072, 4096, 5120, 6144, 8192 ] 
    dt_choice = [ 1e-4, 5e-5, 1e-5, 5e-6, 2e-6, 1e-6, 5e-7, 4e-7, 3e-7, 2e-7, 9e-8, 8e-8, 7e-8, 6e-8, 5e-8, 4e-8, 3e-8, 2e-8, 1e-8 ]
    K_ref = 3.5
    T_width = round(0.1, 2)
    T_step = round(0.05, 2)

    # generate parameter tuples
    tasks = []
    filecount = 0.0
    memcount = 0.0
    for K in K_range:
        for ell1 in ell1_range:
            # for ell2 in ell2_range:
                ell2 = ell1
                elltemp = ell1
                T_target = 0.5*ell1 + 2.0
                # targettemp = 1.02*np.sqrt(1)
                # T_target = 0.5*targettemp + 2.1
                # T_range = np.round(np.array([(T_target / 2) - 3*K ]), 2) # symmetry
                # T_range = np.round(np.array([(2*T_target / 3) - K ]), 2) # init
                T_range = np.round(np.array([-1.43]), 2) # fixed
                # T_range = np.round(np.array([(T_target) - K ]), 2) # Target
                # T_range = np.round(np.array([(T_target - T_width) - K ]), 2) # Lower Bound LB
                # T_range = np.round(np.array([(T_target + T_width) - K ]), 2) # Upper Bound UB
                # T_range = np.round(np.arange((T_target - T_width) - K, (T_target + T_width) - K + T_step/2, T_step), 2) # [LB,UB] branch

                base_index = 2
                tempidx = 11
                K_increment = int(np.floor((K - 4.0) / 0.5 + 1e-10))
                ell_increment = int(np.floor((ell1 - 1.0) / 0.5 + 1e-10))
                choice_index = base_index + K_increment + ell_increment
                N_index = max(0, min(choice_index, len(N_choice) - 1))
                dt_index = max(0, min(choice_index, len(dt_choice) - 1))
                N = N_choice[tempidx] #N_choice[N_index]
                dt = dt_choice[tempidx] #dt_choice[dt_index]
                MPIrank = 128 #max(1, N // 32) # 192 cores per node on Nibi

                for T in T_range:
                    mem_est = 8 * np.ceil((1.1 * np.exp( -4.789714989 + 0.83721882 * np.log10((N**2) * (10.0**T) / dt) - 1.70503490 * elltemp
                                                + 0.23432402 * K + 0.18691756 * elltemp * np.log10((N**2) * (10.0**T) / dt) )) / 8.0 )
                    mem_req = int(np.maximum( 32 , mem_est ))
                    mem = f"{mem_req}G" 
                    addmemcount = mem_req
                    addfilecount = ( N**2 / 4e7 ) * 1e-3 * np.power(10.0,T) / dt
                    memcount += addmemcount
                    filecount += addfilecount
                    tasks.append((float(K), float(ell1), float(ell2), float(T), float(dt), int(N), int(MPIrank), mem))
    return tasks, int(filecount), int(memcount)

# Write runscripts/ and param files
def write_runscripts(tasks, out_dir='runscripts'):
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)

    # timestamp tag to avoid overwriting old param files
    tag = time.strftime("%Y%m%d_%H%M%S")

    # Group tasks by memory and MPI rank count because both are sbatch resources.
    groups = {}
    for idx, t in enumerate(tasks):
        mem = t[-1]
        mpi_ranks = t[-2]
        groups.setdefault((mem, mpi_ranks), []).append((idx, t))

    # write task_params_<mem>_<tag>.txt
    for (mem, mpi_ranks), items in groups.items():
        mem_tag = mem.replace('G', 'G')
        param_fname = out / f"task_params_{mem_tag}_{mpi_ranks}r_{tag}.txt"
        with param_fname.open('w') as pf:
            for (idx, (K, ell1, ell2, T, dt, N, mpi_ranks, mem)) in items:
                pf.write(f"{idx} {K} {ell1:.2f} {ell2:.2f} {T:.2f} "
                         f"{dt:.0e} {N} {mpi_ranks} {mem} "
                         f"{IC} {OPTIMIZE} {TOL} {CONTINUATION} {OPT_T} {SAVE_STATES}\n"
                         )

    # Standalone scripts matching the same C++ invocation used by the array worker.
    for idx, (K, ell1, ell2, T, dt, N, mpi_ranks, mem) in enumerate(tasks):
        mem_gib = int(mem.removesuffix("G"))
        mem_per_cpu_mib = math.ceil(mem_gib * 1024 / mpi_ranks)
        fname = out / f"run_{K}_{ell1:.2f}_{ell2:.2f}_{T:.2f}_{dt:.0e}_{N}_{mpi_ranks}r.sh"
        with fname.open('w') as fh:
            fh.write(f"""#!/bin/bash
#SBATCH --account=def-bprotas
#SBATCH --mail-user=zigicj@mcmaster.ca
#SBATCH --mail-type=ALL
#SBATCH --ntasks={mpi_ranks}
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu={mem_per_cpu_mib}M
#SBATCH --time={SBATCH_TIME}
set -euo pipefail
mkdir -p output
srun ./solver "{IC}" "{N}" "{N}" "{dt:.0e}" "{K}" "{ell1:.2f}" "{ell2:.2f}" "{T:.2f}" \\
    "{OPTIMIZE}" "{TOL}" "{CONTINUATION}" "{OPT_T}" "{SAVE_STATES}" \\
    > "output/run_{idx}_{mpi_ranks}r.log" 2>&1
""")
        os.chmod(fname, 0o755)
    return groups, tag

# Write run_array.sh driver
def write_run_array_sh(groups, tag, out_fname='run_array.sh', max_concurrent=0):
    lines = [
        "#!/bin/bash",
        "set -euo pipefail",
        'mkdir -p slurm_logs output runscripts',
        ''
    ]

    for (mem, mpi_ranks), items in groups.items():
        mem_gib = int(mem.removesuffix("G"))
        mem_per_cpu_mib = math.ceil(mem_gib * 1024 / mpi_ranks)
        param_file = f"runscripts/task_params_{mem}_{mpi_ranks}r_{tag}.txt"
        count = len(items)
        if SEQUENTIAL_TASKS:
            array_spec = f"0-{count-1}%1"  # for sequential tasks
        else:
            array_spec = f"0-{count-1}"  # for concurrent tasks

        lines.append(f"echo \"Group memory={mem}, ranks={mpi_ranks}: tasks={count}, param_file={param_file}\"")

        if SEQUENTIAL_TASKS:
            ### sequential groups ###
            lines.append('if [[ -z "${prev_jobid:-}" ]]; then')
            lines.append(f"    prev_jobid=$(sbatch --parsable --account=def-bprotas --mail-user=zigicj@mcmaster.ca --mail-type=ALL \\")
            lines.append(f"                     --job-name={RUN_ARRAY_NAME} \\")
            lines.append(f"                     --ntasks={mpi_ranks} --cpus-per-task=1 --time={SBATCH_TIME} --mem-per-cpu={mem_per_cpu_mib}M \\")
            lines.append(f"                     --array={array_spec} --output=slurm_logs/slurm-%A_%a.out --error=slurm_logs/slurm-%A_%a.err \\")
            lines.append(f"                     --export=ALL,PARAM_FILE={param_file} ./run_task_array.sh)")
            lines.append("  else")
            lines.append(f"    prev_jobid=$(sbatch --parsable --dependency=afterany:${{prev_jobid}} --account=def-bprotas --mail-user=zigicj@mcmaster.ca --mail-type=ALL \\")
            lines.append(f"                     --job-name={RUN_ARRAY_NAME} \\")
            lines.append(f"                     --ntasks={mpi_ranks} --cpus-per-task=1 --time={SBATCH_TIME} --mem-per-cpu={mem_per_cpu_mib}M \\")
            lines.append(f"                     --array={array_spec} --output=slurm_logs/slurm-%A_%a.out --error=slurm_logs/slurm-%A_%a.err \\")
            lines.append(f"                     --export=ALL,PARAM_FILE={param_file} ./run_task_array.sh)")
            lines.append("  fi")
        else: 
            ### concurrent groups ###
            lines.append(f"prev_jobid=$(sbatch --parsable --account=def-bprotas --mail-user=zigicj@mcmaster.ca --mail-type=ALL \\")
            lines.append(f"                 --job-name={RUN_ARRAY_NAME} \\")
            lines.append(f"                 --ntasks={mpi_ranks} --cpus-per-task=1 --time={SBATCH_TIME} --mem-per-cpu={mem_per_cpu_mib}M \\")
            lines.append(f"                 --array={array_spec} --output=slurm_logs/slurm-%A_%a.out --error=slurm_logs/slurm-%A_%a.err \\")
            lines.append(f"                 --export=ALL,PARAM_FILE={param_file} ./run_task_array.sh)")
        
        lines.append("echo \"Submitted job: $prev_jobid\"")
        lines.append("")

    Path(out_fname).write_text("\n".join(lines))
    os.chmod(out_fname, 0o755)

# Main
def main():
    [tasks,filecount,memcount] = generate_tasks()
    print(f"Parameter combinations (Max 500 concurrent jobs): {len(tasks)}")
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
