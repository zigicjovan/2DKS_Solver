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

# -----------------------------------------------------------------------------
# Select exactly one cluster. Keep this consistent with CLUSTER in the Makefile.
# -----------------------------------------------------------------------------
CLUSTER = "nibi"
# CLUSTER = "anvil"
# CLUSTER = "bridges2"
# CLUSTER = "stampede3"

# Set False to force Anvil wholenode or Bridges-2 RM. When True, Anvil jobs of
# at most 128 ranks use shared and Bridges-2 jobs of at most 64 ranks use
# RM-shared; larger jobs use the corresponding full-node partition.
USE_SHARED_PARTITION = True

CLUSTER_PROFILES = {
    "nibi": {
        "account": "def-bprotas",
        "mail_user": "zigicj@mcmaster.ca",
        "partition": None,
        "shared_partition": None,
        "shared_rank_limit": 0,
        "ranks_per_node": 192,
        "shared_cpus_per_task": 1,
        "shared_memory_gib_per_core": 4,
        "node_memory_gib": 768,
        "time": "02-23:59",
        "modules": (
            "module --force purge",
            "module load StdEnv/2023 gcc/12.3 openmpi/4.1.5 fftw-mpi/3.3.10",
        ),
        "launcher": 'srun --ntasks="{ranks}"',
    },
    "anvil": {
        "account": "mth260100",
        "mail_user": "jovan.zigic@duny.edu",
        "partition": "wholenode",
        "shared_partition": "shared",
        "shared_rank_limit": 128,
        "ranks_per_node": 128,
        "shared_cpus_per_task": 1,
        "shared_memory_gib_per_core": 2,
        "node_memory_gib": 257,
        "time": "02-23:59",
        "modules": (
            "module --force purge",
            "module load gcc/14.2.0 openmpi/4.1.6 fftw/3.3.10",
        ),
        "launcher": 'mpirun --mca pml ucx --mca btl ^openib -np "{ranks}"',
    },
    "bridges2": {
        "account": "mth260024p",
        "mail_user": "jovan.zigic@duny.edu",
        "partition": "RM",
        "shared_partition": "RM-shared",
        "shared_rank_limit": 64,
        "ranks_per_node": 128,
        "shared_cpus_per_task": 1,
        "shared_memory_gib_per_core": 2,
        "node_memory_gib": 256,
        "time": "02-23:59",
        "modules": (
            "module purge",
            "module load gcc/13.3.1-p20240614 openmpi/5.0.8-gcc13.3.1 fftw/3.3.8",
            "hash -r",
        ),
        "launcher": '"$(type -P mpirun)" -np "{ranks}"',
    },
    "stampede3": {
        "account": "tg-mth260100",
        "mail_user": "jovan.zigic@duny.edu",
        "partition": "spr",
        "shared_partition": None,
        "shared_rank_limit": 0,
        "ranks_per_node": 112,
        "shared_cpus_per_task": 1,
        "shared_memory_gib_per_core": None,
        "node_memory_gib": 128,
        "time": "01-23:59",
        "modules": (
            "module reset",
            "module load intel/24.0 impi/21.11 fftw3/3.3.10",
        ),
        "launcher": "ibrun",
    },
}

if CLUSTER not in CLUSTER_PROFILES:
    raise SystemExit(
        f"Unsupported CLUSTER {CLUSTER!r}; use nibi, anvil, bridges2, or stampede3"
    )

PROFILE = CLUSTER_PROFILES[CLUSTER]


def scheduler_settings(mpi_ranks):
    """Return nodes, partition, CPUs/task, and whether to request memory/CPU."""
    shared = (
        USE_SHARED_PARTITION
        and
        PROFILE["shared_partition"] is not None
        and mpi_ranks <= PROFILE["shared_rank_limit"]
    )
    if shared:
        return 1, PROFILE["shared_partition"], PROFILE["shared_cpus_per_task"], False

    nodes = math.ceil(mpi_ranks / PROFILE["ranks_per_node"])
    request_mem_per_cpu = CLUSTER == "nibi"
    return nodes, PROFILE["partition"], 1, request_mem_per_cpu


def allocated_memory_gib(mpi_ranks):
    """Return memory represented by the resources for which the job is charged."""
    nodes, partition, cpus_per_task, _ = scheduler_settings(mpi_ranks)

    # Nibi charges the same CPU cost whether 4 GB per CPU is requested or not.
    if CLUSTER == "nibi":
        return 4 * mpi_ranks

    # Anvil shared and Bridges-2 RM-shared allocate about 2 GB per CPU core.
    if partition == PROFILE["shared_partition"]:
        return math.ceil(
            mpi_ranks * cpus_per_task * PROFILE["shared_memory_gib_per_core"]
        )

    # Anvil wholenode, Bridges-2 RM, and Stampede3 SPR charge whole nodes.
    return nodes * PROFILE["node_memory_gib"]


def memory_gib(mem):
    """Convert values such as '32G' without requiring Python 3.9 removesuffix()."""
    return int(mem[:-1] if mem.endswith("G") else mem)

# User-editable parameter ranges
K_start =           4.0
K_end =             K_start
K_step =            0.5
num_modes_start =   1
num_modes_end =     num_modes_start
ell1_start =        0.00
ell1_end =          0.00
ell1_step =         0.5
ell2_start =        ell1_start
ell2_end =          ell1_end 
ell2_step =         0.1
K_range =           np.round(np.arange(K_start, K_end + K_step/2, K_step), 1)
ell1_range =        10**(np.round(np.arange(ell1_start + 0.01, ell1_end + ell1_step/2, ell1_step), 2))
# ell1_range = [round(x, 2) for x in [np.sqrt(13),np.sqrt(16),np.sqrt(18),np.sqrt(20)]]
ell2_range =        10**(np.round(np.arange(ell2_start + 0.01, ell2_end + ell2_step/2, ell2_step), 2))
# ell2_range = ell1_range

# User-editable global settings
SBATCH_TIME = PROFILE["time"]  # requested time limit (D-HH:MM)
RUN_ARRAY_NAME = f"{CLUSTER}-long-{K_start}-{ell1_start:.2f}_{K_end}-{ell1_end:.2f}.sh"
SEQUENTIAL_TASKS = False

# Solver settings shared by all generated production runs
IC = "s1"
OPTIMIZE = 1
TOL = "1e-10"
CONTINUATION = 1
OPT_T = -3.10
SAVE_STATES = 100

def generate_tasks():
    # Parameter choice formulas for 2DKS problem
    N_choice =  [ 64,   128,  192,  256,  384,  512,  640,  768,  896,  1024, 1536, 2048, 3072, 4096, 5120, 6144, 8192 ] 
    dt_choice = [ 1e-4, 5e-5, 1e-5, 5e-6, 1e-6, 5e-7, 3e-7, 2e-7, 9e-8, 3e-8, 2e-8, 1e-8, 7e-9, 5e-9, 3e-9, 2e-9, 1e-9 ]
    rk_choice = [ 2,    4,    6,    8,    16,   64,   80,   96,   112,  128,  128,  128,  128,  256,  320,  384,  512  ]
    K_ref = 3.5
    T_width = round(0.10, 2)
    T_step = round(0.05, 2)

    # generate parameter tuples
    tasks = []
    filecount = 0.0
    memcount = 0.0
    for K in K_range:
        for ell1 in ell1_range:
            # for ell2 in ell2_range:
                ell2 = ell1
                # elltemp = ell1
                T_target = np.log10( 10**(2.5 - K) * (ell1*ell2)**(0.85) ) 
                # targettemp = 1.02*np.sqrt(1)
                # T_target = 0.5*targettemp + 2.1
                # T_range = np.round(np.array([(2*T_target / 3) - K ]), 2) # init
                # T_range = np.round(np.array([(T_target) - K ]), 2) # Target
                # T_range = np.round(np.array([(T_target - T_width) - K ]), 2) # Lower Bound LB
                # T_range = np.round(np.array([(T_target + T_width) - K ]), 2) # Upper Bound UB
                T_range = np.round(np.arange((T_target - T_width), (T_target + T_width) + T_step/2, T_step), 2) # [LB,UB] branch
                # T_range = np.round(np.array([-2.00]), 2) # fixed

                base_index = 3
                K_increment = int(np.floor((K - 4.0) / 0.5 + 1e-10))
                ell_increment = int(np.floor((ell1 - 1.0) / 2.5 + 1e-10))
                choice_index = base_index + 2*K_increment + ell_increment
                N_index = max(0, min(choice_index, len(N_choice) - 1))
                dt_index = max(0, min(choice_index, len(dt_choice) - 1))
                N = N_choice[N_index-0]
                dt = dt_choice[dt_index]
                MPIrank = rk_choice[N_index-0]
                # tempidx = 8
                # N = 384 #N_choice[tempidx]
                # dt = 7e-9 #dt_choice[tempidx]
                # MPIrank = 128 #max(1, N // 32) # 192 cores per node on Nibi
                nodes, _, _, _ = scheduler_settings(MPIrank)

                for T in T_range:
                    # mem_est = 0.9 * 8 * np.ceil((1.1 * np.exp( -4.789714989 + 0.83721882 * np.log10((N**2) * (10.0**T) / dt) - 1.70503490 * elltemp
                    #                             + 0.23432402 * K + 0.18691756 * elltemp * np.log10((N**2) * (10.0**T) / dt) )) / 8.0 )
                    mem_req = allocated_memory_gib(MPIrank)
                    mem = f"{mem_req}G" 
                    addmemcount = mem_req
                    addfilecount = ( N**2 / 5e7 ) * 1e-3 * np.power(10.0,T) / dt
                    memcount += addmemcount
                    filecount += addfilecount
                    tasks.append((float(K), float(ell1), float(ell2), float(T), float(dt), int(N), int(MPIrank), int(nodes), mem))
    return tasks, int(filecount), int(memcount)

# Write runscripts/ and param files
def write_runscripts(tasks, out_dir='runscripts'):
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)

    # Create a timestamped experiment directory
    tag = time.strftime("%Y%m%d_%H%M%S")
    run_dir = out / tag
    run_dir.mkdir(parents=True, exist_ok=True)

    # Group tasks by memory and MPI rank count because both are sbatch resources.
    groups = {}
    for idx, t in enumerate(tasks):
        mem = t[-1]
        nodes = t[-2]
        mpi_ranks = t[-3]
        groups.setdefault((mem, mpi_ranks, nodes), []).append((idx, t))

    # Write one parameter file per (memory, rank) group
    for (mem, mpi_ranks, nodes), items in groups.items():
        mem_tag = mem.replace('G', 'G')
        param_fname = run_dir / f"task_params_{mem_tag}_{mpi_ranks}r.txt"
        with param_fname.open('w') as pf:
            for (idx, (K, ell1, ell2, T, dt, N, mpi_ranks, nodes, mem)) in items:
                pf.write(f"{idx} {K} {ell1:.2f} {ell2:.2f} {T:.2f} "
                         f"{dt:.0e} {N} {mpi_ranks} {mem} "
                         f"{IC} {OPTIMIZE} {TOL} {CONTINUATION} {OPT_T} {SAVE_STATES}\n"
                         )

    # Standalone scripts matching the same C++ invocation used by the array worker.
    for idx, (K, ell1, ell2, T, dt, N, mpi_ranks, nodes, mem) in enumerate(tasks):
        nodes, partition, cpus_per_task, request_mem_per_cpu = scheduler_settings(mpi_ranks)
        mem_gib = memory_gib(mem)
        mem_per_cpu_mib = math.ceil(mem_gib * 1024 / mpi_ranks)
        ntasks_per_node = min(mpi_ranks, PROFILE["ranks_per_node"])
        partition_line = f"#SBATCH --partition={partition}\n" if partition else ""
        memory_line = f"#SBATCH --mem-per-cpu={mem_per_cpu_mib}M\n" if request_mem_per_cpu else ""
        module_block = "\n".join(PROFILE["modules"])
        launcher = PROFILE["launcher"].format(ranks=mpi_ranks)
        fname = run_dir / f"run_{K}_{ell1:.2f}_{ell2:.2f}_{T:.2f}_{dt:.0e}_{N}_{mpi_ranks}r.sh"
        with fname.open('w') as fh:
            fh.write(f"""#!/bin/bash
#SBATCH --account={PROFILE["account"]}
{partition_line}#SBATCH --job-name={RUN_ARRAY_NAME}
#SBATCH --mail-user={PROFILE["mail_user"]}
#SBATCH --mail-type=ALL
#SBATCH --nodes={nodes}
#SBATCH --ntasks={mpi_ranks}
#SBATCH --ntasks-per-node={ntasks_per_node}
#SBATCH --cpus-per-task={cpus_per_task}
{memory_line}#SBATCH --time={SBATCH_TIME}
set -euo pipefail
{module_block}
mkdir -p output
{launcher} ./solver "{IC}" "{N}" "{N}" "{dt:.0e}" "{K}" "{ell1:.2f}" "{ell2:.2f}" "{T:.2f}" \\
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

    for (mem, mpi_ranks, nodes), items in groups.items():
        nodes, partition, cpus_per_task, request_mem_per_cpu = scheduler_settings(mpi_ranks)
        mem_gib = memory_gib(mem)
        mem_per_cpu_mib = math.ceil(mem_gib * 1024 / mpi_ranks)
        ntasks_per_node = min(mpi_ranks, PROFILE["ranks_per_node"])
        partition_arg = f"--partition={partition}" if partition else ""
        memory_arg = f"--mem-per-cpu={mem_per_cpu_mib}M" if request_mem_per_cpu else ""
        param_file = f"runscripts/{tag}/task_params_{mem}_{mpi_ranks}r.txt"
        count = len(items)
        if SEQUENTIAL_TASKS:
            array_spec = f"0-{count-1}%1"  # for sequential tasks
        elif max_concurrent > 0:
            array_spec = f"0-{count-1}%{max_concurrent}"
        else:
            array_spec = f"0-{count-1}"  # for concurrent tasks

        lines.append(f"echo \"Cluster={CLUSTER}, partition={partition or 'default'}, memory={mem}, ranks={mpi_ranks}: tasks={count}, param_file={param_file}\"")

        if SEQUENTIAL_TASKS:
            ### sequential groups ###
            lines.append('if [[ -z "${prev_jobid:-}" ]]; then')
            lines.append(f"    prev_jobid=$(sbatch --parsable --account={PROFILE['account']} {partition_arg} --mail-user={PROFILE['mail_user']} --mail-type=ALL \\")
            lines.append(f"                     --job-name={RUN_ARRAY_NAME} --nodes={nodes} \\")
            lines.append(f"                     --ntasks={mpi_ranks} --ntasks-per-node={ntasks_per_node} --cpus-per-task={cpus_per_task} --time={SBATCH_TIME} {memory_arg} \\")
            lines.append(f"                     --array={array_spec} --output=slurm_logs/slurm-%A_%a.out --error=slurm_logs/slurm-%A_%a.err \\")
            lines.append(f"                     ./run_task_array.sh {CLUSTER} {param_file})")
            lines.append("else")
            lines.append(f"    prev_jobid=$(sbatch --parsable --dependency=afterany:${{prev_jobid}} --account={PROFILE['account']} {partition_arg} --mail-user={PROFILE['mail_user']} --mail-type=ALL \\")
            lines.append(f"                     --job-name={RUN_ARRAY_NAME} --nodes={nodes} \\")
            lines.append(f"                     --ntasks={mpi_ranks} --ntasks-per-node={ntasks_per_node} --cpus-per-task={cpus_per_task} --time={SBATCH_TIME} {memory_arg} \\")
            lines.append(f"                     --array={array_spec} --output=slurm_logs/slurm-%A_%a.out --error=slurm_logs/slurm-%A_%a.err \\")
            lines.append(f"                     ./run_task_array.sh {CLUSTER} {param_file})")
            lines.append("fi")
        else: 
            ### concurrent groups ###
            lines.append(f"prev_jobid=$(sbatch --parsable --account={PROFILE['account']} {partition_arg} --mail-user={PROFILE['mail_user']} --mail-type=ALL \\")
            lines.append(f"                 --job-name={RUN_ARRAY_NAME} --nodes={nodes} \\")
            lines.append(f"                 --ntasks={mpi_ranks} --ntasks-per-node={ntasks_per_node} --cpus-per-task={cpus_per_task} --time={SBATCH_TIME} {memory_arg} \\")
            lines.append(f"                 --array={array_spec} --output=slurm_logs/slurm-%A_%a.out --error=slurm_logs/slurm-%A_%a.err \\")
            lines.append(f"                 ./run_task_array.sh {CLUSTER} {param_file})")
        
        lines.append("echo \"Submitted job: $prev_jobid\"")
        lines.append("")

    Path(out_fname).write_text("\n".join(lines))
    os.chmod(out_fname, 0o755)

# Main
def main():
    [tasks,filecount,memcount] = generate_tasks()
    print(f"Cluster: {CLUSTER}")
    print(f"Parameter combinations: {len(tasks)}")
    print(f"Filespace required: {(filecount)}K")
    print(f"Memory requested: {(memcount)}GB")
    print(r"File count check: find . -mindepth 2 -type f -printf '%h\n' | sed 's#^\(\./[^/]*\)/.*#\1#' | sort | uniq -c | sort -nr")
    print("Disk usage check: du -h -d 1 | sort -hr")
    groups, tag = write_runscripts(tasks)
    write_run_array_sh(groups, tag, out_fname=RUN_ARRAY_NAME, max_concurrent=0)
    print(f"Experiment directory: runscripts/{tag}")
    print(f"Execution script: bash {RUN_ARRAY_NAME}")
    if CLUSTER == "nibi":
        print("Priority check: sshare -l -A def-bprotas_cpu,rrg-bprotas_cpu -u fu6,zigicj,noahb")
    else:
        print("Queue check: squeue -u $USER")

if __name__ == "__main__":
    main()
