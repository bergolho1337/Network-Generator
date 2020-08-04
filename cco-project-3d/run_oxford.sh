#!/bin/bash

SEEDS=( 1562002891 1562002894 1562005513 1562005553 1562006177 1562007596 1562008172 1562008424 1562009135 1562009769 )

# ========================================================================================================================
# 1) MINIMIZE TERMINALS
# ========================================================================================================================
#PROGRAM_PATH="/home/berg/Github/Network-Generator/cco-project-3d/bin/Cco_3D"
#INPUT_PATH="/home/berg/Github/Network-Generator/cco-project-3d/inputs/05_Lucas/02_Oxford/01_Minimize_Terminals"

#TRANFORM_PROGRAM_PATH="/home/berg/Github/Network-Generator/cco-project-3d/scripts/transform-filter/bin/TransformFilter"
#TRANSFORM_INPUT_PATH="/home/berg/Github/Network-Generator/cco-project-3d/outputs/02_Oxford/01_Minimize_Terminals"
#TRANSFORM_OUTPUT_PATH="/home/berg/Github/MonoAlg3D_C/networks/03_Lucas/02_Oxford/01_Basics/01_Minimize_Terminals"

#for SEED in "${SEEDS[@]}"; do
#    ${PROGRAM_PATH} ${INPUT_PATH}/min\:terminals_seed-${SEED}.ini
#    ${TRANFORM_PROGRAM_PATH} ${TRANSFORM_INPUT_PATH}/seed\:${SEED}/tree_nterm_80.vtk ${TRANSFORM_OUTPUT_PATH}/seed\:${SEED}.vtk 
#done

# ========================================================================================================================
# 2) MINIMIZE TOTAL
# ========================================================================================================================
#PROGRAM_PATH="/home/berg/Github/Network-Generator/cco-project-3d/bin/Cco_3D"
#INPUT_PATH="/home/berg/Github/Network-Generator/cco-project-3d/inputs/05_Lucas/02_Oxford/02_Minimize_Total"

#TRANFORM_PROGRAM_PATH="/home/berg/Github/Network-Generator/cco-project-3d/scripts/transform-filter/bin/TransformFilter"
#TRANSFORM_INPUT_PATH="/home/berg/Github/Network-Generator/cco-project-3d/outputs/02_Oxford/02_Minimize_Total"
#TRANSFORM_OUTPUT_PATH="/home/berg/Github/MonoAlg3D_C/networks/03_Lucas/02_Oxford/01_Basics/02_Minimize_Total"

#for SEED in "${SEEDS[@]}"; do
#    ${PROGRAM_PATH} ${INPUT_PATH}/min\:total_seed-${SEED}.ini
#    ${TRANFORM_PROGRAM_PATH} ${TRANSFORM_INPUT_PATH}/seed\:${SEED}/tree_nterm_80.vtk ${TRANSFORM_OUTPUT_PATH}/seed\:${SEED}.vtk 
#done    

# ========================================================================================================================
# 3) MINIMIZE TERMINALS LINKED
# ========================================================================================================================
PROGRAM_PATH="/home/berg/Github/Network-Generator/cco-project-3d/bin/Cco_3D"
INPUT_PATH="/home/berg/Github/Network-Generator/cco-project-3d/inputs/05_Lucas/02_Oxford/03_Minimize_Terminals_Linked"

TRANFORM_PROGRAM_PATH="/home/berg/Github/Network-Generator/cco-project-3d/scripts/transform-filter/bin/TransformFilter"
TRANSFORM_INPUT_PATH="/home/berg/Github/Network-Generator/cco-project-3d/outputs/02_Oxford/03_Minimize_Terminal_Linked"
TRANSFORM_OUTPUT_PATH="/home/berg/Github/MonoAlg3D_C/networks/03_Lucas/02_Oxford/01_Basics/03_Minimize_Terminal_Linked"

for SEED in "${SEEDS[@]}"; do
    ${PROGRAM_PATH} ${INPUT_PATH}/min\:terminals-linked_seed-${SEED}.ini
    ${TRANFORM_PROGRAM_PATH} ${TRANSFORM_INPUT_PATH}/seed\:${SEED}/tree_nterm_80.vtk ${TRANSFORM_OUTPUT_PATH}/seed\:${SEED}.vtk 
done

# ========================================================================================================================
# 4) MINIMIZE TOTAL LINKED
# ========================================================================================================================
#PROGRAM_PATH="/home/berg/Github/Network-Generator/cco-project-3d/bin/Cco_3D"
#INPUT_PATH="/home/berg/Github/Network-Generator/cco-project-3d/inputs/05_Lucas/02_Oxford/04_Minimize_Total_Linked"

#TRANFORM_PROGRAM_PATH="/home/berg/Github/Network-Generator/cco-project-3d/scripts/transform-filter/bin/TransformFilter"
#TRANSFORM_INPUT_PATH="/home/berg/Github/Network-Generator/cco-project-3d/outputs/02_Oxford/04_Minimize_Total_Linked"
#TRANSFORM_OUTPUT_PATH="/home/berg/Github/MonoAlg3D_C/networks/03_Lucas/02_Oxford/01_Basics/04_Minimize_Total_Linked"

#for SEED in "${SEEDS[@]}"; do
#    ${PROGRAM_PATH} ${INPUT_PATH}/min\:total-linked_seed-${SEED}.ini
#    ${TRANFORM_PROGRAM_PATH} ${TRANSFORM_INPUT_PATH}/seed\:${SEED}/tree_nterm_80.vtk ${TRANSFORM_OUTPUT_PATH}/seed\:${SEED}.vtk 
#done