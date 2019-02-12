#!/bin/bash
PNAME="./bin/LinkedList"

if [ ! -f $PNAME ]; then
	./recompile_project.sh
fi

valgrind --leak-check=full --show-leak-kinds=all ./$PNAME
