#!/bin/bash

set -e  # Stop immediately if any command fails

# Sync solution branches from scratch
rsync -a ~/scratch/Data/SolutionBranches/ \
      ~/projects/rrg-bprotas/zigicj/2DKS_Solver/Data/

# Build aggregate
cd ~/projects/rrg-bprotas/zigicj/2DKS_Solver/
make aggregate

# Generate plots
cd Data/SolutionBranches/
module load gnuplot/6.0.3
gnuplot plot.gp