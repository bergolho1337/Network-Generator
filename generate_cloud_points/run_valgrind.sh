#!/bin/bash
PNAME="./bin/CloudGenerator"
CONFIG_FILE_1="inputs/simple_circle_cloud.ini"
CONFIG_FILE_2="inputs/sparse_circle_cloud.ini"

if [ ! -f $PNAME ]; then
	./recompile_project.sh
fi

valgrind --leak-check=full --show-leak-kinds=all ./$PNAME $CONFIG_FILE_1
