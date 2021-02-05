#!/bin/bash
PNAME="./bin/BifurcationGenerator"

if [ ! -f $PNAME ]; then
	./recompile_project.sh
fi

valgrind --leak-check=full --show-leak-kinds=all ./$PNAME 2 2
