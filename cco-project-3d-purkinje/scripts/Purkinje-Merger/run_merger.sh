#!/bin/bash

SEEDS=( 1562002891 1562002894 1562005513 1562005555 1562006177 )
RATES=( 10 20 40 )

INPUT_PATH="/home/berg/Github/Network-Generator/cco-project-3d-purkinje/scripts/Purkinje-Merger/inputs"
OUTPUT_PATH="/home/berg/Github/Network-Generator/cco-project-3d-purkinje/scripts/Purkinje-Merger/outputs"

for SEED in "${SEEDS[@]}"; do
    for LV_RATE in "${RATES[@]}"; do
		for RV_RATE in "${RATES[@]}"; do
				./bin/Purkinje-Merger inputs/PMJ_Linkrate\:10/His/his_bundle.vtk inputs/PMJ_Linkrate\:${LV_RATE}/LV/LV_seed\:${SEED}.vtk inputs/PMJ_Linkrate\:${RV_RATE}/RV/RV_seed\:${SEED}.vtk outputs/Mixed/seed\:${SEED}_LV_linkrate:${LV_RATE}__RV_linkrate:${RV_RATE}.vtk
		done
    done
    
done
