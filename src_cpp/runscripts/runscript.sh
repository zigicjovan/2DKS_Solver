#!/bin/bash

IC="s1"
N1=64
N2=$N1
dt=1.0e-3
K=1.5
ell1=1.02
ell2=1.02
T=-0.50

optimize=1
tol=1e-6
continuation=1
optT=1
savestates=100

LOGFILE="output/output_IC_${IC}_N1_${N1}_N2_${N2}_dt_${dt}_K_${K}_ell1_${ell1}_ell2_${ell2}.log"

./solver "$IC" "$N1" "$N2" "$dt" "$K" "$ell1" "$ell2" "$T" \
    "$optimize" "$tol" "$continuation" "$optT" "$savestates" 2>&1 | tee "$LOGFILE"