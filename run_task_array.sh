#!/bin/bash
# run_task_array.sh  -- invoked by run_array.sh sbatch
# Uses PARAM_FILE environment variable passed by sbatch --export
# Creates a unique per-job temp dir in /scratch and cleans up at exit.
#
# SBATCH flags (provided by run_array.sh) override the ones below.

#SBATCH --account=def-bprotas
#SBATCH --mail-user=zigicj@mcmaster.ca
#SBATCH --mail-type=ALL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=07-00:00

set -euo pipefail

# SLURM array task id
TASK_ID=${SLURM_ARRAY_TASK_ID:-0}

# PARAM_FILE should be exported by the sbatch command (see run_array.sh)
PARAM_FILE=${PARAM_FILE:-"./runscripts/task_params_32G.txt"}
if [[ ! -f "$PARAM_FILE" ]]; then
    echo "ERROR: parameter file not found: $PARAM_FILE"
    exit 1
fi

# create output and slurm_logs directories
mkdir -p ./output
mkdir -p ./slurm_logs

# create unique scratch temp dir per job
SCRATCH_BASE=${SCRATCH:-/scratch}/${USER}
TMPDIR_JOB="${SCRATCH_BASE}/tmp_job_${SLURM_JOB_ID:-$$}_${SLURM_ARRAY_TASK_ID:-0}"
mkdir -p "$TMPDIR_JOB"
export TMPDIR="$TMPDIR_JOB"
export MATLAB_PREFDIR="$TMPDIR_JOB/matlab_prefs"
mkdir -p "$MATLAB_PREFDIR"
export MFILETMP="$TMPDIR_JOB"

# read the right line (task_id is zero-based; sed is 1-based)
LINE=$(sed -n "$((TASK_ID+1))p" "$PARAM_FILE")
if [[ -z "$LINE" ]]; then
    echo "ERROR: empty line for TASK_ID=$TASK_ID in $PARAM_FILE"
    rm -rf "$TMPDIR_JOB"
    exit 1
fi

# expected line: idx K ell T dt N mem
read IDX K ell T dt N mem <<< "$LINE"

LOG_DIR="./output"
mkdir -p "$LOG_DIR"
ell_str=$(printf "%.2f" "$ell")
T_str=$(printf "%.1f" "$T")
dt_str="$dt"
LOG_FILE="${LOG_DIR}/run_${K}_${ell_str}_${T_str}_${dt_str}_${N}.log"

module load matlab/2024b.1

# MATLAB command for max energy optimization:
MATLAB_CMD="try; main_2DKS(${dt},${N},${K},${K},1,${ell},${ell},0.02,${T},${T},1,'optimize','IC',1e-6,0.0); catch e; disp(getReport(e)); exit(1); end; exit(0);"
# MATLAB command for asymptotic simulations:
# MATLAB_CMD="try; main_2DKS(${dt},${N},${K},${K},1,${ell},${ell},0.02,100,100,1,'plotOptIC','IC',1e-6,${T}); catch e; disp(getReport(e)); exit(1); end; exit(0);"

# --- Retry loop with success-string check ---
attempt=0
max_attempts=2
rc=1

while [[ $attempt -lt $max_attempts ]]; do
    if grep -q "Optimization run complete" "$LOG_FILE" 2>/dev/null; then
        echo "[$(date)] MATLAB reported success. Skipping retry." >> "$LOG_FILE"
        rc=0
        break
    fi

    echo "[$(date)] Running MATLAB attempt $((attempt+1))/$max_attempts for TASK_ID=${TASK_ID} (IDX=${IDX})" >> "$LOG_FILE"
    matlab -nodisplay -nodesktop -nosplash -r "$MATLAB_CMD" >> "$LOG_FILE" 2>&1 || rc=$?

    if grep -q "Optimization run complete" "$LOG_FILE"; then
        echo "[$(date)] MATLAB finished successfully (success string found)." >> "$LOG_FILE"
        rc=0
        break
    else
        echo "[$(date)] MATLAB failed or success string missing. Retrying after 10s..." >> "$LOG_FILE"
        sleep 10
    fi

    attempt=$((attempt+1))
done

# cleanup
rm -rf "$TMPDIR_JOB"

if [[ $rc -ne 0 ]]; then
    echo "[$(date)] TASK ${TASK_ID} finished with error (rc=${rc}). See $LOG_FILE" >&2
    exit $rc
fi

exit 0
