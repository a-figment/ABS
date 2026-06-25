#!/bin/bash
[ ! $# = 1 ] && printf "\033[0;31m$(basename "$0") - Argument error. Require [ STUDY NAME ].\033[0m\n" && exit 1

SOURCE="../Studies/S-${1}/${1}_main.cpp"
OUT="../Studies/S-${1}/run.o"
DBG="../Studies/S-${1}/debug.o"


printf "Bulding '${SOURCE}' ... "
g++ -I ${PWD}/../CppCode/ -Wno-parentheses -Wall -Wextra ${SOURCE} -o ${OUT} -lpthread -Ofast && printf " Done!\n" || exit 1

#g++ -I ${PWD}/../CppCode/ -Wno-parentheses -Wall -Wextra -g ${SOURCE} -o ${DBG} -Ofast && printf " Done!\n" || exit 1

#g++ -I ${PWD}/../CppCode/ -Wno-parentheses -Wall -Wextra ${SOURCE} -g ${OUT} -Ofast 


