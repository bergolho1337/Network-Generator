#!/bin/bash

# Author: Lucas Berg
# This script will generate and transform all the Purkinje networks defined in the INPUT_FOLDER. 
# The output networks will be in the MonoAlg3D domain. 

#SEEDS=( 1562002891 1562002894 1562005513 1562005553 1562006177 1562007596 1562008172 1562008424 1562009134 1562009769 )
SEEDS=( 1562002891 1562002894 1562005513 1562005555 1562006177 )
# ========================================================================================================================
# 1) MINIMIZE LENGTH (linking PMJ's at the end)
# ========================================================================================================================
PROGRAM_PATH="/home/berg/Github/Network-Generator/cco-project-3d-purkinje/bin/Cco_3D"
INPUT_PATH="/home/berg/Github/Network-Generator/cco-project-3d-purkinje/inputs/03_PMJ_Linktries:20"
OUTPUT_PATH="/home/berg/Github/Network-Generator/cco-project-3d-purkinje/outputs/03_PMJ_Linktries:20"
PURKINJE_MERGER_PATH="/home/berg/Github/Network-Generator/cco-project-3d-purkinje/scripts/Purkinje-Merger"
PURKINJE_MERGER_INPUT_PATH="/home/berg/Github/Network-Generator/cco-project-3d-purkinje/scripts/Purkinje-Merger/inputs/04_PMJ_Linkrate:40"
PURKINJE_MERGER_OUTPUT_PATH="/home/berg/Github/Network-Generator/cco-project-3d-purkinje/scripts/Purkinje-Merger/outputs/04_PMJ_Linkrate:40"
ELECTRIC_ERROR_PATH="/home/berg/Github/Network-Generator/cco-project-3d-purkinje/scripts/Electric-Calculator"
ELECTRIC_ERROR_INPUT_PATH="/home/berg/Github/Network-Generator/cco-project-3d-purkinje/scripts/Electric-Calculator/inputs/PMJ_Linkrate:40"
ELECTRIC_ERROR_OUTPUT_PATH="/home/berg/Github/Network-Generator/cco-project-3d-purkinje/scripts/Electric-Calculator/outputs/PMJ_Linkrate:40"
#MONOALG3D_PATH="/home/berg/Github/MonoAlg3D_C"

# LEFT VENTRICLE
#for SEED in "${SEEDS[@]}"; do
#    ${PROGRAM_PATH} ${INPUT_PATH}/elizabeth_biventricular_coupled_co:length_seed:${SEED}_LV.ini & 
#done

# RIGHT VENTRICLE
#for SEED in "${SEEDS[@]}"; do
#    ${PROGRAM_PATH} ${INPUT_PATH}/elizabeth_biventricular_coupled_co:length_seed:${SEED}_RV_back_top.ini ${INPUT_PATH}/elizabeth_biventricular_coupled_co:length_seed:${SEED}_RV_front_top.ini ${INPUT_PATH}/elizabeth_biventricular_coupled_co:length_seed:${SEED}_RV_front_bottom.ini ${INPUT_PATH}/elizabeth_biventricular_coupled_co:length_seed:${SEED}_RV_back_bottom.ini ${OUTPUT_PATH}/RV_seed:${SEED}
#done

# Copying
#for SEED in "${SEEDS[@]}"; do
#    cp ${OUTPUT_PATH}/LV_seed:${SEED}/tree_nterm_650.vtk ${PURKINJE_MERGER_INPUT_PATH}/LV/LV_seed:${SEED}.vtk
#    cp ${OUTPUT_PATH}/RV_seed:${SEED}/tree_nterm_192.vtk ${PURKINJE_MERGER_INPUT_PATH}/RV/RV_seed:${SEED}.vtk
#done

# Merging
#for SEED in "${SEEDS[@]}"; do
#    ${PURKINJE_MERGER_PATH}/bin/Purkinje-Merger ${PURKINJE_MERGER_INPUT_PATH}/His/his_bundle.vtk ${PURKINJE_MERGER_INPUT_PATH}/LV/LV_seed:${SEED}.vtk ${PURKINJE_MERGER_INPUT_PATH}/RV/RV_seed:${SEED}.vtk ${PURKINJE_MERGER_OUTPUT_PATH}/LVRV_seed:${SEED}.vtk
#done

# Electric error
#for SEED in "${SEEDS[@]}"; do
#    ${ELECTRIC_ERROR_PATH}/bin/PurkinjeError ${ELECTRIC_ERROR_PATH}/inputs/elizabeth_gold_standart_LVRV.vtk ${ELECTRIC_ERROR_INPUT_PATH}/LVRV_seed:${SEED}.vtk ${ELECTRIC_ERROR_PATH}/inputs/elizabeth_pmjs_LVRV.txt ${ELECTRIC_ERROR_OUTPUT_PATH}/LVRV_LAT_seed:${SEED}.vtk
#done

# ========================================================================================================================
# 2) MINIMIZE LENGTH (linking PMJ's at intervals)
# ========================================================================================================================
