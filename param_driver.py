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

# ----------------------------
# User-editable parameter ranges
# ----------------------------
K_range = range(0, 6)
ell_range = np.round(np.arange(1.06, 1.51, 0.08), 2)

# ----------------------------
# Generate parameter tuples
# ----------------------------
def generate_tasks():
    tasks = []
    for K in K_range:
        for ell in ell_range:
            if K == 0:
                T_range = np.round(np.arange(0.0, 2.1, 0.1), 1)
                dt, N, mem = 1e-3, 24, '32G'
            elif K == 1:
                T_range = np.round(np.arange(-0.5, 1.6, 0.1), 1)
                dt, N, mem = 1e-3, 32, '32G'
            elif K == 2:
                T_range = np.round(np.arange(-1.0, 1.1, 0.1), 1)
                dt, N, mem = 1e-3, 32, '32G'
            elif K == 3:
                T_range = np.round(np.arange(-1.5, 0.6, 0.1), 1)
                dt, N, mem = 5e-5, 48, '32G'
            elif K == 4:
                T_range = np.round(np.arange(-2.0, -0.9, 0.1), 1)
                dt, N, mem = 1e-6, 64, '32G'
            elif K == 5:
                T_range = np.round(np.arange(-3.0, -1.9, 0.1), 1)
                dt, N, mem = 5e-8, 128, '64G'
            else:
                continue

            for T in T_range:
                tasks.append((K, float(ell), float(T), dt, int(N), mem))
    return tasks

# ----------------------------
# Write runscripts/ and param files
# ----------------------------
def write_runscripts(tasks, out_dir='runscripts'):
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)

    # timestamp tag to avoid overwriting old param files
    tag = time.strftime("%Y%m%d_%H%M%S")

    # group tasks by memory
    groups = {}
    for idx, t in enumerate(tasks):
        mem = t[5]
        groups.setdefault(mem, []).append((idx, t))

    # write task_params_<mem>_<tag>.txt
    for mem, items in groups.items():
        mem_tag = mem.replace('G', 'G')
        param_fname = out / f"task_params_{mem_tag}_{tag}.txt"
        with param_fname.open('w') as pf:
            for (idx, (K, ell, T, dt, N, mem)) in items:
                output_file = f"output/run_{K}_{ell:.2f}_{T:.1f}_{dt}_{N}.mat"
                pf.write(f"{idx} {K} {ell:.2f} {T:.1f} {dt:.0e} {N} {mem} {output_file}\n")

    # placeholder scripts
    for idx, (K, ell, T, dt, N, mem) in enumerate(tasks):
        fname = out / f"run_{K}_{ell:.2f}_{T:.1f}_{dt}_{N}.sh"
        with fname.open('w') as fh:
            fh.write(f"""#!/bin/bash
#SBATCH --account=def-bprotas
#SBATCH --mail-user=zigicj@mcmaster.ca
#SBATCH --mail-type=ALL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem={mem}
#SBATCH --time=07-00:00
module load matlab/2024b.1
# Placeholder; real execution uses ../run_task_array.sh
""")
    return groups, tag

# ----------------------------
# Write run_array.sh driver
# ----------------------------
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
        lines.append(f"         --ntasks=1 --cpus-per-task=8 --time=07-00:00 --mem={mem} \\")
        lines.append(f"         --array={array_spec} --output=slurm_logs/slurm-%A_%a.out --error=slurm_logs/slurm-%A_%a.err \\")
        lines.append(f"         --export=ALL,PARAM_FILE={param_file} ./run_task_array.sh")
        lines.append("else")
        lines.append(f"  echo \"Dry run: would sbatch --array={array_spec} --mem={mem} --cpus-per-task=8 --export=PARAM_FILE={param_file} ./run_task_array.sh\"")
        lines.append("fi\n")

    Path(out_fname).write_text("\n".join(lines))
    os.chmod(out_fname, 0o755)

# ----------------------------
# Main
# ----------------------------
def main():
    tasks = generate_tasks()
    print(f"Total parameter combinations: {len(tasks)}")
    groups, tag = write_runscripts(tasks)
    write_run_array_sh(groups, tag, max_concurrent=0)
    print(f"Wrote runscripts/ (task param files with tag {tag}) and run_array.sh")

if __name__ == "__main__":
    main()
