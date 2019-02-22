#!/bin/bash
PNAME="./bin/Cco"
QPERF="1"
PPERF="10"
PTERM="1"
RPERF="1"
NTERM="10"

if [ ! -f $PNAME ]; then
	./recompile_project.sh
fi

valgrind --leak-check=full --show-leak-kinds=all ./$PNAME $QPERF $PPERF $PTERM $RPERF $NTERM
