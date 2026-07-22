#!/bin/bash

IC="s1"
N1=64
N2=$N1
dt=1.0e-4
K=0.5
ell1=1.02
ell2=$ell1
T=-1.50

optimize=1
tol=1e-6
continuation=1
optT=1
savestates=100
MPI_RANKS=1

LOGFILE="output/output_IC_${IC}_N1_${N1}_N2_${N2}_dt_${dt}_K_${K}_ell1_${ell1}_ell2_${ell2}_opt_${optimize}_cont_${continuation}_optT_${optT}_rank_${MPI_RANKS}.log"

mpirun -np "$MPI_RANKS" ./solver "$IC" "$N1" "$N2" "$dt" "$K" "$ell1" "$ell2" "$T" \
    "$optimize" "$tol" "$continuation" "$optT" "$savestates" 2>&1 | tee "$LOGFILE"
