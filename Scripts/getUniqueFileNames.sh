#!/bin/bash

# Keys are uniquely identified by "--" character
FILE_LIST=("$@")
FILE_NAMES=()
for f in ${FILE_LIST[@]}; do
    if [[ $f == *"--"* ]]; then
        TEMP=${f#*--}  # Remove prefix
        FILE_NAME=${TEMP%.*}  # Remove suffix 
        FILE_NAMES+=(${FILE_NAME})
    fi
done

UNIQUE_FILE_NAMES=($(printf '%s\n' "${FILE_NAMES[@]}" | sort -u))
printf '%s\n' "${UNIQUE_FILE_NAMES[@]}"
#printf "No. of files = ${#UNIQUE_FILE_NAMES[@]}. \n"

