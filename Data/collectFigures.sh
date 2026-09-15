#!/usr/bin/env bash
set -euo pipefail

# Run from the directory containing the testcase directories.
# Override with PYTHON=/path/to/venv/bin/python if needed.
script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
python_bin="${PYTHON:-python3}"
shopt -s nullglob
testcases=()

for testcase_dir in */; do
    testcase="${testcase_dir%/}"
    files=("${testcase}"/ForwardSolution/fwd_*.dat)
    energy_files=("${testcase}"/EnergyEvolution/energy*.dat)
    gif="forwardSolution/movie${testcase}.gif"

    if [[ -f "$gif" ]]; then
        echo "Skipping ${testcase}: ${gif} already exists."
        continue
    fi

    if ((${#energy_files[@]} == 0)); then
        echo "Skipping ${testcase}: no EnergyEvolution output yet."
        continue
    fi

    if ((${#files[@]})); then
        testcases+=("$testcase")
    fi
done

if ((${#testcases[@]} == 0)); then
    echo 'No unfinished testcase directories found.' >&2
    exit 0
fi

exec "$python_bin" "$script_dir/postprocess/generateFigures.py" "${testcases[@]}" "$@"
