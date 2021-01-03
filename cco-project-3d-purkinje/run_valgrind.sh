#!/bin/bash
PNAME="./bin/Cco_3D"
CONFIG_FILE="inputs/test_LV.ini"
#CONFIG_FILE="inputs/test_RV_back_top.ini inputs/test_RV_front_top.ini inputs/test_RV_front_bottom.ini inputs/test_RV_back_bottom.ini"

if [ ! -f $PNAME ]; then
	./recompile_project.sh
fi

valgrind --leak-check=full --show-leak-kinds=all $PNAME $CONFIG_FILE
