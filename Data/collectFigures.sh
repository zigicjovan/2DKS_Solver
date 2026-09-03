#!/usr/bin/env bash
set -euo pipefail

# Run from the directory containing the testcase directories.
# Override with PYTHON=/path/to/venv/bin/python if needed.
script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
python_bin="${PYTHON:-python3}"
shopt -s nullglob
testcases=()
for testcase_dir in */; do
    files=("${testcase_dir%/}"/ForwardSolution/fwd_*.dat)
    if ((${#files[@]})); then
        testcases+=("${testcase_dir%/}")
    fi
done
if ((${#testcases[@]} == 0)); then
    echo 'No testcase directories containing ForwardSolution/fwd_*.dat found.' >&2
    exit 1
fi
exec "$python_bin" "$script_dir/postprocess/generateFigures.py" "${testcases[@]}" "$@"
