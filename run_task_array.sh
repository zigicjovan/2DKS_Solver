#!/bin/bash
# run_task_array.sh  -- invoked by run_array.sh sbatch
# Receives the cluster name and parameter file as positional arguments.
# Launches one C++/MPI solver run for the selected parameter row.

set -euo pipefail

# Manual default used only when this script is invoked directly. The generated
# array driver passes its selected cluster explicitly.
DEFAULT_CLUSTER="nibi"
# DEFAULT_CLUSTER="anvil"
# DEFAULT_CLUSTER="bridges2"
# DEFAULT_CLUSTER="stampede3"

CLUSTER=${1:-$DEFAULT_CLUSTER}

# SLURM array task id
TASK_ID=${SLURM_ARRAY_TASK_ID:-0}

# --- Stagger start of each array task (3 sec per task, modulo 60 sec) ---
sleep $(( (TASK_ID * 3) % 60 ))

# Passing this as an argument avoids sbatch --export portability problems.
PARAM_FILE=${2:-"./runscripts/task_params_32G.txt"}
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

# expected line: idx K ell1 ell2 T dt N MPI_ranks mem IC optimize tol continuation optT savestates
read -r IDX K ell1 ell2 T dt N MPI_RANKS mem \
    IC_str optimize tol continuation optT savestates <<< "$LINE"

case "$CLUSTER" in
    nibi)
        module --force purge
        module load StdEnv/2023 gcc/12.3 openmpi/4.1.5 fftw-mpi/3.3.10
        LAUNCH=(srun --ntasks="$MPI_RANKS")
        ;;
    anvil)
        module --force purge
        module load gcc/14.2.0 openmpi/4.1.6 fftw/3.3.10
        LAUNCH=(mpirun --mca pml ucx --mca btl '^openib' -np "$MPI_RANKS")
        ;;
    bridges2)
        module purge
        module load gcc/13.3.1-p20240614
        module load openmpi/5.0.8-gcc13.3.1
        module load fftw/3.3.8
        hash -r

        if [[ -z "${LOCAL:-}" ]]; then
            echo "ERROR: Bridges-2 LOCAL directory is unavailable." >&2
            exit 1
        fi

        export SLURM_TMPDIR="$LOCAL"

        MPI_LAUNCHER=$(type -P mpirun)
        if [[ -z "$MPI_LAUNCHER" ]]; then
            echo "ERROR: Bridges-2 mpirun executable not found." >&2
            exit 1
        fi

        LAUNCH=("$MPI_LAUNCHER" -np "$MPI_RANKS")
        ;;
    stampede3)
        module reset
        module load intel/24.0
        module load impi/21.11
        module load fftw3/3.3.10

        export SLURM_TMPDIR="/tmp/$USER/job_${SLURM_JOB_ID}"
        mkdir -p "$SLURM_TMPDIR"

        LAUNCH=(ibrun)
        ;;
    *)
        echo "ERROR: unsupported cluster '$CLUSTER'." >&2
        exit 2
        ;;
esac

LOG_DIR="./output"
mkdir -p "$LOG_DIR"
ell1_str=$(printf "%.2f" "$ell1")
ell2_str=$(printf "%.2f" "$ell2")
T_str=$(printf "%.2f" "$T")
dt_str="$dt"
LOG_FILE="${LOG_DIR}/run_${IC_str}_${K}_${ell1_str}_${ell2_str}_${T_str}_${dt_str}_${N}_${MPI_RANKS}r.log"

# --- Write SLURM Job ID to log file ---
echo -e "\n=============================" >> "$LOG_FILE"
echo "SLURM_JOB_ID: ${SLURM_JOB_ID:-N/A}" >> "$LOG_FILE"
echo "SLURM_ARRAY_JOB_ID: ${SLURM_ARRAY_JOB_ID:-N/A}  SLURM_ARRAY_TASK_ID: ${SLURM_ARRAY_TASK_ID:-N/A}" >> "$LOG_FILE"
echo "CLUSTER: $CLUSTER" >> "$LOG_FILE"
echo "=============================" >> "$LOG_FILE"

# --- Retry once if the solver returns a nonzero exit code ---
attempt=0
max_attempts=2
rc=1

while [[ $attempt -lt $max_attempts ]]; do
    echo "[$(date)] Running C++ attempt $((attempt+1))/$max_attempts for TASK_ID=${TASK_ID} (IDX=${IDX})" >> "$LOG_FILE"
    rc=0
    "${LAUNCH[@]}" ./solver "$IC_str" "$N" "$N" "$dt" "$K" "$ell1" "$ell2" "$T" \
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
