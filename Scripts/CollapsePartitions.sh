#!/bin/bash

[ $# = 0 ] && printf "\033[0;31m$(basename "$0") - Argument error. Require [ STUDY NAME, NO. PARTITIONS ].\033[0m\n" && exit 1


# 1. Get all file names in the study
# 2. 
