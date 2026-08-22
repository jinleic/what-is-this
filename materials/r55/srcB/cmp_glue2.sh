#!/bin/bash
# Compare glue2 full-stratum output multiset against the already-produced outB files.
set -e
cd /Users/jinleic/jinleic-workspace/math/r55/srcB
for spec in "13 52" "16 71"; do
  n=${spec% *}; e=${spec#* }
  rm -f /tmp/g2_all.g6 /tmp/g2_all.csv
  dmin=$(( (2*e + n - 1) / n ))
  dmax=$(( n-1 < 13 ? n-1 : 13 ))
  for d in $(seq "$dmin" "$dmax"); do
    q=$((n-1-d)); [ "$q" -lt 0 ] && continue
    ./glue2 "$n" "$e" "$d" 0 9999999 /tmp/g2_$d.g6 /tmp/g2_$d.csv ../data 2>/dev/null
    cat /tmp/g2_$d.g6 >> /tmp/g2_all.g6
    cat /tmp/g2_$d.csv >> /tmp/g2_all.csv
  done
  echo "($n,$e):"
  if diff <(sort /tmp/g2_all.g6) <(sort ../outB/r45_${n}_${e}.g6) >/dev/null; then
    echo "  g6 multiset identical"
  else echo "  G6 DIFF"; fi
  if diff <(sort /tmp/g2_all.csv) <(tail -n +2 ../outB/r45_${n}_${e}.counts.csv | sort) >/dev/null; then
    echo "  counts identical"
  else echo "  COUNTS DIFF"; fi
done
