#!/bin/bash

SEEDS=( 1562002891 1562002894 1562005513 1562005555 1562006177 )
RATES=( 10 20 40 )

INPUT_PATH="/home/berg/Github/Network-Generator/cco-project-3d-purkinje/scripts/Electric-Calculator/inputs"
OUTPUT_PATH="/home/berg/Github/Network-Generator/cco-project-3d-purkinje/scripts/Electric-Calculator/outputs"

for SEED in "${SEEDS[@]}"; do
    for LV_RATE in "${RATES[@]}"; do
		for RV_RATE in "${RATES[@]}"; do
			./bin/PurkinjeError inputs/elizabeth_gold_standart_LVRV.vtk inputs/Mixed/seed\:${SEED}_LV_linkrate\:${LV_RATE}__RV_linkrate\:${RV_RATE}.vtk inputs/elizabeth_pmjs_LVRV.txt outputs/Mixed/seed\:${SEED}_LV_linkrate:${LV_RATE}__RV_linkrate:${RV_RATE}.vtk
        done
    done
done
