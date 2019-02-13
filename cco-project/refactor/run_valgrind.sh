#!/bin/bash
PNAME="./bin/Cco"

if [ ! -f $PNAME ]; then
	./recompile_project.sh
fi

valgrind --leak-check=full --show-leak-kinds=all ./$PNAME
