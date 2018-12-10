#!/bin/bash
PNAME="./bin/FDMMonodomain1D"

if [ "$#" -ne 1 ]; then
	echo "[ERROR] Illegal number of parameters"
	exit 1
fi

if [ ! -f $PNAME ]; then
	./recompile_project.sh
fi

CONFIG_FILEPATH=$1

valgrind --leak-check=full --show-leak-kinds=all ./$PNAME $CONFIG_FILEPATH
