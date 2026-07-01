#!/bin/bash

IC="s1"
N1=128
N2=$N1
dt=1e-4
K=3
ell1=1.02
ell2=1.02
T=-1.00

optimize=1
tol=1e-6
continuation=0
optT=-3

./solver $IC $N1 $N2 $dt $K $ell1 $ell2 $T $optimize $tol $continuation $optT