#!/bin/bash

set -e  # Stop immediately if any command fails

# Sync solution branches from scratch
rsync -a $SCRATCH/Data/SolutionBranches/ \
      $PROJECT/x-jzigic/2DKS_Solver/Data/

# Build aggregate
cd $PROJECT/x-jzigic/2DKS_Solver/
make aggregate

# Generate plots
cd Data/SolutionBranches/
module load gnuplot
gnuplot plot.gp
