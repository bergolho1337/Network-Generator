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
INPUT_PATH="/home/berg/Github/Network-Generator/cco-project-3d-purkinje/inputs/01_PMJ_Final"
OUTPUT_PATH="/home/berg/Github/Network-Generator/cco-project-3d-purkinje/outputs/01_PMJ_Final"
PURKINJE_MERGER_PATH="/home/berg/Github/Network-Generator/cco-project-3d-purkinje/scripts/Purkinje-Merger"
#MONOALG3D_PATH="/home/berg/Github/MonoAlg3D_C"

# LEFT VENTRICLE
#for SEED in "${SEEDS[@]}"; do
#    ${PROGRAM_PATH} ${INPUT_PATH}/elizabeth_biventricular_coupled_co:length_seed:${SEED}_LV.ini & 
#done

# RIGHT VENTRICLE
for SEED in "${SEEDS[@]}"; do
    ${PROGRAM_PATH} ${INPUT_PATH}/elizabeth_biventricular_coupled_co:length_seed:${SEED}_RV_back_top.ini ${INPUT_PATH}/elizabeth_biventricular_coupled_co:length_seed:${SEED}_RV_front_top.ini ${INPUT_PATH}/elizabeth_biventricular_coupled_co:length_seed:${SEED}_RV_front_bottom.ini ${INPUT_PATH}/elizabeth_biventricular_coupled_co:length_seed:${SEED}_RV_back_bottom.ini ${OUTPUT_PATH}/RV_seed:${SEED}
    #cp ${OUTPUT_PATH}/full_RV/tree_nterm_192.vtk ${OUTPUT_PATH}/RV/elizabeth_total_length_seed:${SEED}_RV.vtk
done


# LEFT VENTRICLE (copying)
#for SEED in "${SEEDS[@]}"; do
#    cp ${OUTPUT_PATH}/LV/elizabeth_total_length_seed:${SEED}_LV/tree_nterm_650.vtk ${OUTPUT_PATH}/LV/elizabeth_total_length_seed:${SEED}_LV.vtk
#done

# RIGHT VENTRICLE (linking)
#for SEED in "${SEEDS[@]}"; do
#    cp ${OUTPUT_PATH}/RV/elizabeth_total_length_seed:${SEED}_RV_back_bottom/tree_nterm_50.vtk ${TERMINAL_LINKER_PATH}/inputs/elizabeth_min_length_back_bottom.vtk
#    cp ${OUTPUT_PATH}/RV/elizabeth_total_length_seed:${SEED}_RV_back_top/tree_nterm_10.vtk ${TERMINAL_LINKER_PATH}/inputs/elizabeth_min_length_back_top.vtk
#    cp ${OUTPUT_PATH}/RV/elizabeth_total_length_seed:${SEED}_RV_front_bottom/tree_nterm_150.vtk ${TERMINAL_LINKER_PATH}/inputs/elizabeth_min_length_front_bottom.vtk
#    cp ${OUTPUT_PATH}/RV/elizabeth_total_length_seed:${SEED}_RV_front_top/tree_nterm_5.vtk ${TERMINAL_LINKER_PATH}/inputs/elizabeth_min_length_front_top.vtk
#    ${TERMINAL_LINKER_PATH}/bin/TerminalLinker ${TERMINAL_LINKER_PATH}/inputs/elizabeth_min_length_back_top.vtk ${TERMINAL_LINKER_PATH}/inputs/elizabeth_min_length_back_bottom.vtk ${TERMINAL_LINKER_PATH}/inputs/elizabeth_min_length_front_top.vtk ${TERMINAL_LINKER_PATH}/inputs/elizabeth_min_length_front_bottom.vtk ${OUTPUT_PATH}/RV/elizabeth_total_length_seed:${SEED}_RV.vtk
#done

# LVRV (merging)
#for SEED in "${SEEDS[@]}"; do
#    cp ${OUTPUT_PATH}/LV/elizabeth_total_length_seed:${SEED}_LV.vtk ${PURKINJE_MERGER_PATH}/inputs/01_CO_Length/LV/elizabeth_total_length_seed:${SEED}_LV.vtk
#    cp ${OUTPUT_PATH}/RV/elizabeth_total_length_seed:${SEED}_RV.vtk ${PURKINJE_MERGER_PATH}/inputs/01_CO_Length/RV/elizabeth_total_length_seed:${SEED}_RV.vtk
#    ${PURKINJE_MERGER_PATH}/bin/Purkinje-Merger ${PURKINJE_MERGER_PATH}/inputs/01_CO_Length/His/elizabeth_min_length_his.vtk ${PURKINJE_MERGER_PATH}/inputs/01_CO_Length/LV/elizabeth_total_length_seed:${SEED}_LV.vtk ${PURKINJE_MERGER_PATH}/inputs/01_CO_Length/RV/elizabeth_total_length_seed:${SEED}_RV.vtk ${PURKINJE_MERGER_PATH}/outputs/01_CO_Length/elizabeth_total_length_seed:${SEED}_LVRV.vtk
#    cp ${PURKINJE_MERGER_PATH}/outputs/01_CO_Length/elizabeth_total_length_seed:${SEED}_LVRV.vtk ${MONOALG3D_PATH}/networks/elizabeth_total_length_seed:${SEED}_LVRV.vtk
#done

# ========================================================================================================================
# 2) MINIMIZE LENGTH (linking PMJ's at intervals)
# ========================================================================================================================
