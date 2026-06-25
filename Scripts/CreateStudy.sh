#!/bin/bash
[ ! $# = 1 ] && printf "\033[0;31m$(basename "$0") - Argument error. Require [ STUDY NAME ].\033[0m\n" && exit 1

#while getopts ":s:" arg; do
#  case $arg in
#    s) STUDY_NAME=$OPTARG;;
#  esac
#done

STUDY_NAME=$1
STUDY_DIR="S-${STUDY_NAME}"

[[ -d ../Studies/$STUDY_DIR ]] && printf "\033[0;31m$(basename "$0") - Study already exists! Exiting ... \033[0m\n" && exit 1

mkdir -p ../Studies/${STUDY_DIR}/Data
cp ../Studies/S-template/template_main.cpp ../Studies/${STUDY_DIR}/${STUDY_NAME}_main.cpp

touch ../Studies/${STUDY_DIR}/Study-Description.txt && echo "An optional 'reminder' file for information on what the Study is about." >> ../Studies/${STUDY_DIR}/Study-Description.md

printf "\033[0;32mStudy \"$STUDY_NAME\" created!\033[0m\n"

