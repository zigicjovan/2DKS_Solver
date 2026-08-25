#!/usr/bin/env bash
set -euo pipefail

project_data_dir="$HOME/projects/rrg-bprotas/zigicj/2DKS_Solver/Data"
scratch_data_dir="$HOME/scratch/Data"

cd "$project_data_dir"

testcase_dirs=()

for testcase_dir in "$scratch_data_dir"/*/; do
    [[ -d "$testcase_dir" ]] || continue

    forward_dir="${testcase_dir%/}/ForwardSolution"

    if compgen -G "${forward_dir}/f*.dat" > /dev/null; then
        testcase_dirs+=("$testcase_dir")
    fi
done

total=${#testcase_dirs[@]}
echo "Found ${total} testcase director$([[ "$total" -eq 1 ]] && echo "y" || echo "ies") to process."

for index in "${!testcase_dirs[@]}"; do
    testcase_dir="${testcase_dirs[$index]}"
    forward_dir="${testcase_dir%/}/ForwardSolution"
    testcase="$(basename "${testcase_dir%/}")"

    current=$((index + 1))

    echo
    echo "${current} of ${total}: Generating figures for ${testcase}"

    matlab -batch \
        "set(0,'DefaultFigureVisible','off'); generateFigures('${testcase}')"

    find "$forward_dir" -maxdepth 1 -type f -name 'f*.dat' -delete

    echo "${current} of ${total}: Deleted forward data."
done
