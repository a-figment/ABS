#!/bin/bash
[ $# = 0 ] && printf "\033[0;31m$(basename "$0") - Argument error. Require [ STUDY NAME, DATA KEY INDICES (OPT) ].\033[0m\n" && exit 1

STUDY_NAME=$1
STUDY_DIR="S-${STUDY_NAME}"
SOURCE="${STUDY_NAME}_main.cpp"
SOURCE_LOC="../Scripts/${STUDY_DIR}/"
DATA_DIR=../Studies/${STUDY_DIR}/Data

# Optional re-build of source code
read -p "Build study? [Enter] " key
if [[ $key = "" ]]; then
    source Build.sh $STUDY_NAME #"${SOURCE_LOC}/${SOURCE}"
fi

read -p "Clean study? [Enter] " key
if [[ $key = "" ]]; then
    rm ${DATA_DIR}/*dat &>/dev/null
    printf "\033[0;32m\nData files removed from Study $STUDY_NAME.\033[0m\n" 
fi

#./CleanStudy.sh ${STUDY_NAME} - need to allow for repeat runs

printf "\033[0;32mRunning study \"$STUDY_NAME\" ... \033[0m\n"
ARGS=("$@")

cd ../Studies/$STUDY_DIR/ && ./run.o ${ARGS[@]}

