#!/bin/bash

# -----------------------------
# Dry-run option: execute "bash run_array.sh --dry-run" and tune param_driver.py to adjust
# -----------------------------
DRY_RUN=0
MAX_CONCURRENT=0   # 0 = no concurrency limit
if [[ "$1" == "--dry-run" ]]; then
    DRY_RUN=1
    echo "Dry run mode: no jobs will be submitted."
fi

# -----------------------------
# Load modules
# -----------------------------
module load matlab/2024b.1
module load python/3.11

# -----------------------------
# Helper functions
# -----------------------------
count_params() {
    group=$1
    python3 -c "import param_driver; print(param_driver.total_params($group))"
}

get_mem() {
    group=$1
    python3 -c "import param_driver; params=param_driver.get_params($group); print(params[0][-1] if params else '0')"
}

available_groups() {
    python3 -c "import param_driver; print(' '.join(map(str,param_driver.available_array_groups())))"
}

# -----------------------------
# Submit function
# -----------------------------
submit_array() {
    group=$1
    total=$(count_params $group)
    mem=$(get_mem $group)

    if [[ $total -eq 0 ]]; then
        echo "Array group $group: no tasks to submit"
        return
    fi

    # Decide array specification
    if [[ $MAX_CONCURRENT -gt 0 ]]; then
        ARRAY_SPEC="0-$(($total - 1))%$MAX_CONCURRENT"
    else
        ARRAY_SPEC="0-$(($total - 1))"
    fi

    # Print example MATLAB command for first task
    if [[ "$DRY_RUN" -eq 1 ]]; then
        python3 param_driver.py --array-group=$group --dry-run
    fi

    if [[ "$DRY_RUN" -eq 0 ]]; then
        sbatch --account=def-bprotas \
               --mail-user=zigicj@mcmaster.ca \
               --mail-type=ALL \
               --ntasks=1 \
               --cpus-per-task=8 \
               --time=07-00:00 \
               --mem=$mem \
               --array=$ARRAY_SPEC <<EOF
#!/bin/bash
module load matlab/2024b.1
module load python/3.11
python3 param_driver.py \$SLURM_ARRAY_TASK_ID --array-group=$group
EOF
    else
        echo "Dry run: would submit array group $group with --array=$ARRAY_SPEC --mem=$mem"
    fi
}

# -----------------------------
# Submit all available array groups
# -----------------------------
for group in $(available_groups); do
    submit_array $group
done
