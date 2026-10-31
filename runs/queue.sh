#!/bin/bash
cd "$(dirname "$0")/.."
export OMP_NUM_THREADS=1
while IFS= read -r line; do
  [ -z "$line" ] && continue
  echo "$line | $(python3 membrane_fem.py $line 2>&1 | tail -1)" >> "runs/$2"
done < "runs/$1"
echo DONE >> "runs/$2"
