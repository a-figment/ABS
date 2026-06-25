#!/bin/bash
[ ! $# = 1 ] && printf "\033[0;31m$(basename "$0") - Argument error. Require [ STUDY NAME ].\033[0m\n" && exit 1

STUDY_NAME=$1
STUDY_DIR="S-${STUDY_NAME}"
DATA_DIR=../Studies/${STUDY_DIR}/Data

rm ${DATA_DIR}/*dat &>/dev/null
printf "\033[0;32m\nData files removed from Study $STUDY_NAME.\033[0m\n" 


#FTR=$(ls $DATA_DIR/*dat &>/dev/null)
#FTR=$(ls ${DATA_DIR}/*dat)

#printf "\033[0;32m\nData files removed from Study $STUDY_NAME:\033[0m\n" 
#echo $FTR | tr ' ' '\n' | sort -V | head -n5
#echo -e " .\n .\n . "
#echo $FTR | tr ' ' '\n' | sort -V | tail -n5
#
#printf "\n${FTR[@]}"
#printf "\n${#FTR[@]}\n"
#
#if [ ${#FTR[@]} -gt 0 ]; then
#    rm ${DATA_DIR}/*dat 
#    printf "\033[0;32m\nData files removed from Study $STUDY_NAME:\033[0m\n" 
#    echo $FTR | tr ' ' '\n' | sort -V | head -n5
#    echo -e " .\n .\n . "
#    echo $FTR | tr ' ' '\n' | sort -V | tail -n5
#else
#    printf "\033[0;32mNo data files to remove from Study $STUDY_NAME.\033[0m\n"
#fi


#[[ ! -f "$DATA_DIR/*dat" ]] && printf "\033[0;32mNo data files to remove from Study $STUDY_NAME.\033[0m\n"


