#!/bin/bash

[[ $# = 1 ]] && OUT=$1 || OUT="run.o"

# Run 
g++ -Wno-parentheses -Wall -Wextra _main.cpp -o  $OUT -lpthread -Ofast

# Debug
g++ -g _main.cpp -o debug.o -lpthread

