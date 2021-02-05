#!/bin/bash

#SEEDS=( 1562002891 1562002894 1562005513 1562005553 1562008172 )	# Patient-Specific (Experiment 1)
#SEEDS=( 1562002894 1562005513 1562005553 1562006177 1562008172 )	# Patient-Specific (Experiment 2)
SEEDS=( 1562002891 1562002894 1562005513 1562007596 1562008172 )	# Patient-Specific (Experiment 3)

INPUT_PATH="/home/berg/Github/Network-Generator/cco-project-3d-purkinje/scripts/Purkinje-Merger/inputs"
OUTPUT_PATH="/home/berg/Github/Network-Generator/cco-project-3d-purkinje/scripts/Purkinje-Merger/outputs"

for SEED in "${SEEDS[@]}"; do
	#./bin/Purkinje-Merger inputs/02_Patient_Specific/His/His.vtk inputs/02_Patient_Specific/Experiment-1/LV_seed:${SEED}.vtk inputs/02_Patient_Specific/Experiment-1/RV_seed:${SEED}.vtk outputs/02_Patient_Specific/Experiment-1/LVRV_seed:${SEED}.vtk
	#./bin/Purkinje-Merger inputs/02_Patient_Specific/His/His.vtk inputs/02_Patient_Specific/Experiment-2/LV_seed:${SEED}.vtk inputs/02_Patient_Specific/Experiment-2/RV_seed:${SEED}.vtk outputs/02_Patient_Specific/Experiment-2/LVRV_seed:${SEED}.vtk
	./bin/Purkinje-Merger inputs/02_Patient_Specific/His/His.vtk inputs/02_Patient_Specific/Experiment-3/LV_seed:${SEED}.vtk inputs/02_Patient_Specific/Experiment-3/RV_seed:${SEED}.vtk outputs/02_Patient_Specific/Experiment-3/LVRV_seed:${SEED}.vtk
done
