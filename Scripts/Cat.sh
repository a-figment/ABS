#!/bin/bash

#[ $# = 0 ] && printf "\033[0;31m$(basename "$0") - Argument error. Require [ STUDY NAME, DATA KEY INDICES (OPT) ].\033[0m\n" && exit 1


# Takes a file name, and some number of key elements. All of the missing key elements are globbed and concatenated to a file
# e.g. ./Cat.sh MFPT s0 b1 - Concatenates the files s0-b1-p*--MFPT.dat. This case is trivial since key elements are in order.

ARGS=("$@")
FILE_NAME=${ARGS[0]}
DROP_KEYS=${ARGS[1]}

FILES=$(ls *${DROP_KEYS}*${FILE_NAME}*)
N_KEYS=$((${#ARGS[@]} - 1))
KEY=""

for i in $(seq 2 $N_KEYS); do
    FILES=$(ls ${FILES} | grep "${ARGS[${i}]}")
    KEY+="${ARGS[${i}]}-"
done
KEY+="-"

printf "Concatenating ... "
printf "Done! New file \033[0;34m${KEY}${FILE_NAME}.dat\033[0m produced. Key element dropped: ${DROP_KEYS}. \n"
#printf '%s\n' ${FILES[@]}

cat ${FILES[@]} > "${KEY}${FILE_NAME}.dat"

#printf "Done! New file ${KEY}${FILE_NAME}.dat produced. \n"

