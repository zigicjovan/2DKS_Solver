#!/bin/bash

mkdir -p ./SolutionBranches/EnergyEvolution

for dir in */; do
    dir=${dir%/}

    if [ -f "$dir/EnergyEvolution/energy.dat" ]; then
        cp "$dir/EnergyEvolution/energy.dat" \
           "./SolutionBranches/EnergyEvolution/${dir}.dat"
    fi
done
