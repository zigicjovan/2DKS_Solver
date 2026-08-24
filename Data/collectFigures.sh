#!/usr/bin/env bash
set -euo pipefail

# project_data_dir="$HOME/projects/rrg-bprotas/zigicj/2DKS_Solver/Data"
# scratch_data_dir="$HOME/scratch/Data"

project_data_dir="$HOME/Desktop/2DKS/2DKS_cpp/Data"
scratch_data_dir="$HOME/Desktop/2DKS/2DKS_cpp/Data"

cd "$project_data_dir"

for testcase_dir in "$scratch_data_dir"/*/; do
    [[ -d "$testcase_dir" ]] || continue

    forward_dir="${testcase_dir%/}/ForwardSolution"

    # Skip cases with no forward solution files.
    if ! compgen -G "${forward_dir}/f*.dat" > /dev/null; then
        continue
    fi

    testcase="$(basename "${testcase_dir%/}")"

    echo "Generating figures for: ${testcase}"

    # Assumes your MATLAB function is generateFigures(testcase).
    matlab -batch "set(0,'DefaultFigureVisible','off'); generateFigures('${testcase}')"

    # Runs only if MATLAB exited successfully.
    find "$forward_dir" -maxdepth 1 -type f -name 'f*.dat' -delete

    echo "Deleted forward data for: ${testcase}"
done
