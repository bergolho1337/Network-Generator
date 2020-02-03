#!/bin/bash
PNAME="./bin/DLA"

if [ ! -f $PNAME ]; then
	./recompile_project.sh
fi

valgrind --leak-check=full --show-leak-kinds=all ./$PNAME inputs/square_walker.ini 
#valgrind ./$PNAME inputs/square_walker.ini 
