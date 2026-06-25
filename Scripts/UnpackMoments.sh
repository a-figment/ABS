#!/bin/bash

NROWS=$1; shift
OUT=$1; shift
FILES=( "$@" )

echo "Files to process: "
#printf "%s\n" "${FILES[@]}"
echo "%s\n" "${FILES[@]}"
echo "Number of moments = ${NROWS}\n"
echo "Output file = ${OUT}_x.dat\n"

# Every row is a new file.
for i in `seq 1 1 ${NROWS}`; do
    for j in "${FILES[@]}"; do
        echo "${j}"
        tail -n+${i} ${j} | head -1 >> "${OUT}_${i}.dat"
    done
done
