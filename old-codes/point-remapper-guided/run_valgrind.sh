#!/bin/bash
PNAME="./bin/graph"

if [ "$#" -ne 1 ]; then
	echo "[ERROR] Illegal number of parameters"
	exit 1
fi

if [ ! -f $PNAME ]; then
	./recompile_project.sh
fi

INPUT_FILEPATH=$1

#valgrind --leak-check=full --show-leak-kinds=all --track-origins=no ./$PNAME $INPUT_FILEPATH
valgrind --leak-check=full --show-leak-kinds=all --undef-value-errors=no ./$PNAME $INPUT_FILEPATH
