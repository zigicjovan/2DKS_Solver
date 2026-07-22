#!/bin/bash
# run_task_array.sh  -- invoked by run_array.sh sbatch
# Uses PARAM_FILE environment variable passed by sbatch --export
# Launches one C++/MPI solver run for the selected parameter row.

set -euo pipefail

# SLURM array task id
TASK_ID=${SLURM_ARRAY_TASK_ID:-0}

# --- Stagger start of each array task (3 sec per task, modulo 60 sec) ---
sleep $(( (TASK_ID * 3) % 60 ))

# PARAM_FILE should be exported by the sbatch command (see run_array.sh)
PARAM_FILE=${PARAM_FILE:-"./runscripts/task_params_32G.txt"}
if [[ ! -f "$PARAM_FILE" ]]; then
    echo "ERROR: parameter file not found: $PARAM_FILE"
    exit 1
fi

# create output and slurm_logs directories
mkdir -p ./output
mkdir -p ./slurm_logs

# read the right line (task_id is zero-based; sed is 1-based)
LINE=$(sed -n "$((TASK_ID+1))p" "$PARAM_FILE")
if [[ -z "$LINE" ]]; then
    echo "ERROR: empty line for TASK_ID=$TASK_ID in $PARAM_FILE"
    exit 1
fi

# expected line: idx K ell1 ell2 T dt N MPI_ranks mem
read -r IDX K ell1 ell2 T dt N MPI_RANKS mem \
    IC_str optimize tol continuation optT savestates <<< "$LINE"

LOG_DIR="./output"
mkdir -p "$LOG_DIR"
ell1_str=$(printf "%.2f" "$ell1")
ell2_str=$(printf "%.2f" "$ell2")
T_str=$(printf "%.2f" "$T")
dt_str="$dt"
LOG_FILE="${LOG_DIR}/run_${IDX}_${IC_str}_${K}_${ell1_str}_${ell2_str}_${T_str}_${dt_str}_${N}_${MPI_RANKS}r.log"

# --- Write SLURM Job ID to log file ---
echo -e "\n=============================" >> "$LOG_FILE"
echo "SLURM_JOB_ID: ${SLURM_JOB_ID:-N/A}" >> "$LOG_FILE"
echo "SLURM_ARRAY_JOB_ID: ${SLURM_ARRAY_JOB_ID:-N/A}  SLURM_ARRAY_TASK_ID: ${SLURM_ARRAY_TASK_ID:-N/A}" >> "$LOG_FILE"
echo "=============================" >> "$LOG_FILE"

# --- Retry once if the solver returns a nonzero exit code ---
attempt=0
max_attempts=2
rc=1

module load fftw-mpi/3.3.10

while [[ $attempt -lt $max_attempts ]]; do
    echo "[$(date)] Running C++ attempt $((attempt+1))/$max_attempts for TASK_ID=${TASK_ID} (IDX=${IDX})" >> "$LOG_FILE"
    rc=0
    srun --ntasks="$MPI_RANKS" ./solver "$IC_str" "$N" "$N" "$dt" "$K" "$ell1" "$ell2" "$T" \
        "$optimize" "$tol" "$continuation" "$optT" "$savestates" >> "$LOG_FILE" 2>&1 || rc=$?

    if [[ $rc -eq 0 ]]; then
        echo "[$(date)] C++ solver finished successfully." >> "$LOG_FILE"
        rc=0
        break
    else
        echo "[$(date)] C++ solver failed with rc=${rc}." >> "$LOG_FILE"
        sleep 10
    fi

    attempt=$((attempt+1))
done

if [[ $rc -ne 0 ]]; then
    echo "[$(date)] TASK ${TASK_ID} finished with error (rc=${rc}). See $LOG_FILE" >&2
    exit $rc
fi
