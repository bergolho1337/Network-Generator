#!/bin/bash

# Author: Lucas Berg
# This script will generate and transform all the Purkinje networks defined in the INPUT_FOLDER. 
# The output networks will be in the MonoAlg3D domain. 

SEEDS=( 1562002891 1562002894 1562005513 1562005553 1562006177 1562007596 1562008172 1562008424 1562009134 1562009769 )

# ========================================================================================================================
# 1) MINIMIZE LENGTH
# ========================================================================================================================
PROGRAM_PATH="/home/berg/Github/Network-Generator/cco-project-3d/bin/Cco_3D"
INPUT_PATH="/home/berg/Github/Network-Generator/cco-project-3d/inputs/05_Lucas/04_Elizabeth/01_CO_Length"

TRANFORM_PROGRAM_PATH="/home/berg/Github/Network-Generator/cco-project-3d/scripts/transform-filter/bin/TransformFilter"
TRANSFORM_INPUT_PATH="/home/berg/Github/Network-Generator/cco-project-3d/outputs/04_Elizabeth/01_CO_Length"
TRANSFORM_OUTPUT_PATH="/home/berg/Documentos/Purkinje-Merger/inputs/01_CO_Length"

#for SEED in "${SEEDS[@]}"; do
#    ${PROGRAM_PATH} ${INPUT_PATH}/elizabeth_biventricular_coupled_co:length_seed:${SEED}_RV.ini & 
#done

#for SEED in "${SEEDS[@]}"; do
#    ${PROGRAM_PATH} ${INPUT_PATH}/elizabeth_biventricular_coupled_co:length_seed:${SEED}_LV.ini & 
#done

#for SEED in "${SEEDS[@]}"; do
#    ${TRANFORM_PROGRAM_PATH} ${TRANSFORM_INPUT_PATH}/RV_seed\:${SEED}/tree_nterm_140.vtk ${TRANSFORM_OUTPUT_PATH}/RV/RV_seed\:${SEED}.vtk
#done

#for SEED in "${SEEDS[@]}"; do
#    ${TRANFORM_PROGRAM_PATH} ${TRANSFORM_INPUT_PATH}/LV_seed\:${SEED}/tree_nterm_650.vtk ${TRANSFORM_OUTPUT_PATH}/LV/LV_seed\:${SEED}.vtk
#done

# ========================================================================================================================
# 2) MINIMIZE TOTAL ACTIVATION TIME
# ========================================================================================================================
PROGRAM_PATH="/home/berg/Github/Network-Generator/cco-project-3d/bin/Cco_3D"
INPUT_PATH="/home/berg/Github/Network-Generator/cco-project-3d/inputs/05_Lucas/04_Elizabeth/02_CO_Activation_Time"

TRANFORM_PROGRAM_PATH="/home/berg/Github/Network-Generator/cco-project-3d/scripts/transform-filter/bin/TransformFilter"
TRANSFORM_INPUT_PATH="/home/berg/Github/Network-Generator/cco-project-3d/outputs/04_Elizabeth/02_CO_Activation_Time"
TRANSFORM_OUTPUT_PATH="/home/berg/Documentos/Purkinje-Merger/inputs/02_CO_Activation_Time"

#for SEED in "${SEEDS[@]}"; do
#    ${PROGRAM_PATH} ${INPUT_PATH}/elizabeth_biventricular_coupled_co:at_seed:${SEED}_RV.ini & 
#done

#for SEED in "${SEEDS[@]}"; do
#    ${PROGRAM_PATH} ${INPUT_PATH}/elizabeth_biventricular_coupled_co:at_seed:${SEED}_LV.ini & 
#done

#for SEED in "${SEEDS[@]}"; do
#    ${TRANFORM_PROGRAM_PATH} ${TRANSFORM_INPUT_PATH}/RV_seed\:${SEED}/tree_nterm_140.vtk ${TRANSFORM_OUTPUT_PATH}/RV/RV_seed\:${SEED}.vtk
#done

#for SEED in "${SEEDS[@]}"; do
#    ${TRANFORM_PROGRAM_PATH} ${TRANSFORM_INPUT_PATH}/LV_seed\:${SEED}/tree_nterm_650.vtk ${TRANSFORM_OUTPUT_PATH}/LV/LV_seed\:${SEED}.vtk
#done

# ========================================================================================================================
# 3) MINIMIZE TERMINAL ACTIVATION TIME
# ========================================================================================================================
PROGRAM_PATH="/home/berg/Github/Network-Generator/cco-project-3d/bin/Cco_3D"
INPUT_PATH="/home/berg/Github/Network-Generator/cco-project-3d/inputs/05_Lucas/04_Elizabeth/03_CO_Terminals"

TRANFORM_PROGRAM_PATH="/home/berg/Github/Network-Generator/cco-project-3d/scripts/transform-filter/bin/TransformFilter"
TRANSFORM_INPUT_PATH="/home/berg/Github/Network-Generator/cco-project-3d/outputs/04_Elizabeth/03_CO_Terminals"
TRANSFORM_OUTPUT_PATH="/home/berg/Documentos/Purkinje-Merger/inputs/03_CO_Terminals"

#for SEED in "${SEEDS[@]}"; do
#    ${PROGRAM_PATH} ${INPUT_PATH}/elizabeth_biventricular_coupled_co:at_seed:${SEED}_RV.ini & 
#done

for SEED in "${SEEDS[@]}"; do
    ${PROGRAM_PATH} ${INPUT_PATH}/elizabeth_biventricular_coupled_co:at_seed:${SEED}_LV.ini & 
done

#for SEED in "${SEEDS[@]}"; do
#    ${TRANFORM_PROGRAM_PATH} ${TRANSFORM_INPUT_PATH}/RV_seed\:${SEED}/tree_nterm_140.vtk ${TRANSFORM_OUTPUT_PATH}/RV/RV_seed\:${SEED}.vtk
#done

#for SEED in "${SEEDS[@]}"; do
#    ${TRANFORM_PROGRAM_PATH} ${TRANSFORM_INPUT_PATH}/LV_seed\:${SEED}/tree_nterm_650.vtk ${TRANSFORM_OUTPUT_PATH}/LV/LV_seed\:${SEED}.vtk
#done
