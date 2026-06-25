#!/bin/bash
[ $# = 0 ] && printf "\033[0;31m$(basename "$0") - Argument error. Require [ STUDY NAME, DATA KEY INDICES (OPT) ].\033[0m\n" && exit 1


STUDY_NAME=$1
STUDY_DIR="S-${STUDY_NAME}"
SOURCE="${STUDY_NAME}_main.cpp"
SOURCE_LOC="../Scripts/${STUDY_DIR}/"

# Optional re-build of source code
printf "\033[0;32mBuilding study \"$STUDY_NAME\" ... \033[0m\n"
source Build.sh $STUDY_NAME #"${SOURCE_LOC}/${SOURCE}"


printf "\033[0;32mCleaning study \"$STUDY_NAME\" ... \033[0m\n"
./CleanStudy.sh ${STUDY_NAME}

printf "\033[0;32mRunning study \"$STUDY_NAME\" ... \033[0m\n"
ARGS=("$@")

cd ../Studies/$STUDY_DIR/ && ./run.o ${ARGS[@]}

