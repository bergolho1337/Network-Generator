#!/bin/bash
PNAME="./bin/Cco"
CONFIG_FILE="inputs/simple_cco.ini"

if [ ! -f $PNAME ]; then
	./recompile_project.sh
fi

valgrind --leak-check=full --show-leak-kinds=all ./$PNAME $CONFIG_FILE
