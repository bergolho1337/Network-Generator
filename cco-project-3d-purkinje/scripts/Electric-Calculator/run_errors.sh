#!/bin/bash

#SEEDS=( 1562002891 1562002894 1562005513 1562005554 1562006177 1562007596 1562008172 1562008424 1562009134 1562009769 ) 
SEEDS=( 1562002891 1562002894 1562005513 1562005554 1562006177 )

PROGRAM_PATH="/home/berg/Documentos/Elizabeth/Scripts/Purkinje-Error"
#INPUT_FOLDER="/home/berg/Documentos/Elizabeth/Scripts/Purkinje-Error/inputs/01_CO_Length"
#INPUT_FOLDER="/home/berg/Documentos/Elizabeth/Scripts/Purkinje-Error/inputs/02_CO_Length_PMJ_last"
#INPUT_FOLDER="/home/berg/Documentos/Elizabeth/Scripts/Purkinje-Error/inputs/03_CO_Length_Multiobjective"
INPUT_FOLDER="/home/berg/Documentos/Elizabeth/Scripts/Purkinje-Error/inputs/04_CO_Length_PMJ_interval_with_region_radius"

for SEED in "${SEEDS[@]}"; do
    echo "=========================================================================================================================="
    ${PROGRAM_PATH}/bin/PurkinjeError inputs/elizabeth_gold_standart_LVRV.vtk ${INPUT_FOLDER}/elizabeth_total_length_seed\:${SEED}_LVRV.vtk inputs/elizabeth_pmjs_LVRV.txt
    echo "=========================================================================================================================="
done
