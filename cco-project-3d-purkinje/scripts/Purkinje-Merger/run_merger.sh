#!/bin/bash

SEEDS=( 1562002891 1562002894 1562005513 1562005553 1562006177 1562007596 1562008172 1562008424 1562009134 1562009769 ) 

INPUT_PATH="/home/berg/Documentos/Purkinje-Merger/inputs/01_CO_Length"
OUTPUT_PATH="/home/berg/Documentos/Purkinje-Merger/outputs/01_CO_Length"

for SEED in "${SEEDS[@]}"; do
    ./bin/Purkinje-Merger ${INPUT_PATH}/His/His.vtk ${INPUT_PATH}/LV/LV_seed:${SEED}.vtk ${INPUT_PATH}/RV/RV_seed:${SEED}.vtk ${OUTPUT_PATH}/LVRV_seed:${SEED}.vtk 
done

INPUT_PATH="/home/berg/Documentos/Purkinje-Merger/inputs/02_CO_Activation_Time"
OUTPUT_PATH="/home/berg/Documentos/Purkinje-Merger/outputs/02_CO_Activation_Time"

for SEED in "${SEEDS[@]}"; do
    ./bin/Purkinje-Merger ${INPUT_PATH}/His/His.vtk ${INPUT_PATH}/LV/LV_seed:${SEED}.vtk ${INPUT_PATH}/RV/RV_seed:${SEED}.vtk ${OUTPUT_PATH}/LVRV_seed:${SEED}.vtk 
done
