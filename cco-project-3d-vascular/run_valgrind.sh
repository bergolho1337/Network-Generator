#!/bin/bash
PNAME="./bin/Cco_3D"
CONFIG_FILE="inputs/01_Basics/minimize_volume_using_cloud_points.ini"

if [ ! -f $PNAME ]; then
	./recompile_project.sh
fi

valgrind --leak-check=full --show-leak-kinds=all $PNAME $CONFIG_FILE
