#!/bin/bash
[ $# = 0 ] && printf "\033[0;31m$(basename "$0") - Argument error. Require [ STUDY NAME, DATA KEY INDICES (OPT) ].\033[0m\n" && exit 1

STUDY_NAME=$1
JOB_ID=$2
PARTITION_ID=$3

STUDY_DIR="S-${STUDY_NAME}"
SOURCE="${STUDY_NAME}_main.cpp"

cd "../Studies/${STUDY_DIR}" 
ARGS=("$@")

printf "\033[0;32mRunning study \"$STUDY_NAME\" ... \033[0m\n"
./run.o ${ARGS[@]}



