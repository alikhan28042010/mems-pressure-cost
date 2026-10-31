#!/bin/bash
cd "$(dirname "$0")"
mkdir -p runs
for args in \
  "--fine 8 --coarse 30" \
  "--fine 5 --coarse 25" \
  "--fine 4 --coarse 25" \
  "--fine 8 --coarse 30 --frame 150 100" \
  "--fine 6 --coarse 30 --frame 150 100" \
  "--fine 5 --coarse 30 --frame 150 100" \
  "--fine 6 --coarse 30 --frame 150 100 --t 9" \
  "--fine 6 --coarse 30 --frame 150 100 --t 11" \
  "--fine 6 --coarse 30 --frame 150 100 --t 15" \
  "--fine 6 --coarse 30 --frame 150 100 --t 17"; do
  name=$(echo $args | tr -d ' -')
  python3 ccx_membrane.py $args --out runs/$name
done
