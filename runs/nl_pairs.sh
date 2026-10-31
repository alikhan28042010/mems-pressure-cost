#!/bin/bash
cd "$(dirname "$0")/.."
export OMP_NUM_THREADS=1
out=runs/resNL.txt
s="--p 50 --pi 3 -1 80 --resL 60 10 10 --resT 60 10 20 --nxy 8 --nz 1"
for pair in "12 1097" "9 1093"; do
  set -- $pair
  echo "LIN --t $1 --side $2 $s | $(python3 membrane_fem.py --t $1 --side $2 $s 2>&1 | tail -1)" >> $out
  echo "NL --t $1 --side $2 $s | $(python3 nonlinear_check.py --t $1 --side $2 $s --steps 4 2>&1 | tail -1)" >> $out
done
echo "LIN su13 | $(python3 membrane_fem.py --nxy 8 --nz 1 2>&1 | tail -1)" >> $out
echo "NL su13 | $(python3 nonlinear_check.py --nxy 8 --nz 1 --steps 1 2>&1 | tail -1)" >> $out
echo DONE >> $out
