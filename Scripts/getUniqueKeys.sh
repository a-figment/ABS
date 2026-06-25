#!/bin/bash
# Keys are uniquely identified by "--" character
FILE_LIST=("$@")
KEYS=()
for f in ${FILE_LIST[@]}; do
    if [[ $f == *"--"* ]]; then
        TEMP=${f#*--}  # Remove prefix
        FILE_NAME=${TEMP%.*}  # Remove suffix 
        KEYS+=(${f%--*})
    fi
done
UNIQUE_KEYS=($(printf '%s\n' "${KEYS[@]}" | sort -u))
printf '%s\n' "${UNIQUE_KEYS[@]}"
#printf "No. of unique keys = ${#UNIQUE_KEYS[@]}. \n"


