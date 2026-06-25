#!/bin/bash

[ ! $# = 2 ] && printf "\033[0;31m$(basename "$0") - Argument error. Require [ FILE NAME, NUM INSTANCES ].\033[0m\n" && exit 1

#
FILE_ID=$1

#
N_VERSIONS=$2

echo $PWD

printf "Concatenating $N_VERSIONS instances of $FILE_ID ... "
for i in $(seq 0 ${N_VERSIONS}); do
    cat p*-${FILE_ID}_${i}.dat > ${FILE_ID}_${i}.dat
done
printf "Done!\n"

DIR="pX-${FILE_ID}--Raw/"
read -p "Move raw data to $DIR? [y/n] " key
[[ $key = "y" ]] && mv p*-${FILE_ID}_*.dat $DIR
