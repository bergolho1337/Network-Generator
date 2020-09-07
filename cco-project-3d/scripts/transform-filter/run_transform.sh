# ===============================================================================================================
# Apply a pre-defined Transform filter using VTK over a set of Purkinje Network that are inside a folder
# =============================================================================================================== 

#!/bin/bash

SEEDS=( 1562002891 1562002894 1562005513 1562005553 1562006177 1562007596 1562008172 1562008424 1562009134 1562009769 ) 
INPUT_PATH="/home/berg/Github/Network-Generator/cco-project-3d/outputs/04_Elizabeth/02_CO_Activation_Time"
OUTPUT_PATH="/home/berg/Documentos/Purkinje-Merger/inputs/02_CO_Activation_Time/RV"

for SEED in "${SEEDS[@]}"; do
    ./bin/TransformFilter ${INPUT_PATH}/RV_seed:${SEED}/tree_nterm_140.vtk ${OUTPUT_PATH}/RV_seed:${SEED}.vtk 
done

