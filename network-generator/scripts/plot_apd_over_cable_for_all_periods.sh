#!/bin/bash

START_PERIOD=195
END_PERIOD=100

BCLS=($(seq $START_PERIOD -5 $END_PERIOD))

for BCL in "${BCLS[@]}"; do
	FOLDER_NAME=$BCL"ms"
	echo "[!] Plotting folder $BCL ..."
	python calc_apd_di_directory.py $FOLDER_NAME $BCL
done
