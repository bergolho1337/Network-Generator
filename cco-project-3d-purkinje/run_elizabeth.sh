#!/bin/bash

# Author: Lucas Berg
# This script will generate and transform all the Purkinje networks defined in the INPUT_FOLDER. 
# The output networks will be in the MonoAlg3D domain. 

#SEEDS=( 1562002891 1562002894 1562005513 1562005553 1562006177 1562007596 1562008172 1562008424 1562009134 1562009769 )
SEEDS=( 1562002891 1562002894 1562005513 1562005554 1562006177 )
# ========================================================================================================================
# 1) MINIMIZE LENGTH (default)
# ========================================================================================================================
#PROGRAM_PATH="/home/berg/Documentos/Elizabeth/Scripts/Network-Generator/bin/Cco_3D"
#INPUT_PATH="/home/berg/Documentos/Elizabeth/Scripts/Network-Generator/inputs/01_CO_Length"
#OUTPUT_PATH="/home/berg/Documentos/Elizabeth/Scripts/Network-Generator/outputs/01_CO_Length"
#TERMINAL_LINKER_PATH="/home/berg/Documentos/Elizabeth/Scripts/Terminal-Linker"
#PURKINJE_MERGER_PATH="/home/berg/Documentos/Elizabeth/Scripts/Purkinje-Merger"
#MONOALG3D_PATH="/home/berg/Github/MonoAlg3D_C"

# LEFT VENTRICLE
#for SEED in "${SEEDS[@]}"; do
#    ${PROGRAM_PATH} ${INPUT_PATH}/elizabeth_biventricular_coupled_co:length_seed:${SEED}_LV.ini & 
#done

# RIGHT VENTRICLE (back bottom)
#for SEED in "${SEEDS[@]}"; do
#    ${PROGRAM_PATH} ${INPUT_PATH}/elizabeth_biventricular_coupled_co:length_seed:${SEED}_RV_back_bottom.ini & 
#done

# RIGHT VENTRICLE (back top)
#for SEED in "${SEEDS[@]}"; do
#    ${PROGRAM_PATH} ${INPUT_PATH}/elizabeth_biventricular_coupled_co:length_seed:${SEED}_RV_back_top.ini & 
#done

# RIGHT VENTRICLE (front bottom)
#for SEED in "${SEEDS[@]}"; do
#    ${PROGRAM_PATH} ${INPUT_PATH}/elizabeth_biventricular_coupled_co:length_seed:${SEED}_RV_front_bottom.ini & 
#done

# RIGHT VENTRICLE (front top)
#for SEED in "${SEEDS[@]}"; do
#    ${PROGRAM_PATH} ${INPUT_PATH}/elizabeth_biventricular_coupled_co:length_seed:${SEED}_RV_front_top.ini & 
#done

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
# 2) MINIMIZE LENGTH (linking PMJ's at the end)
# ========================================================================================================================
#PROGRAM_PATH="/home/berg/Documentos/Elizabeth/Scripts/Network-Generator/bin/Cco_3D"
#INPUT_PATH="/home/berg/Documentos/Elizabeth/Scripts/Network-Generator/inputs/02_CO_Length_PMJ_last"
#OUTPUT_PATH="/home/berg/Documentos/Elizabeth/Scripts/Network-Generator/outputs/02_CO_Length_PMJ_last"
#TERMINAL_LINKER_PATH="/home/berg/Documentos/Elizabeth/Scripts/Terminal-Linker"
#PURKINJE_MERGER_PATH="/home/berg/Documentos/Elizabeth/Scripts/Purkinje-Merger"
#MONOALG3D_PATH="/home/berg/Github/MonoAlg3D_C"

# LEFT VENTRICLE
#for SEED in "${SEEDS[@]}"; do
#    ${PROGRAM_PATH} ${INPUT_PATH}/elizabeth_biventricular_coupled_co:length_seed:${SEED}_LV.ini & 
#done

# RIGHT VENTRICLE (back bottom)
#for SEED in "${SEEDS[@]}"; do
#    ${PROGRAM_PATH} ${INPUT_PATH}/elizabeth_biventricular_coupled_co:length_seed:${SEED}_RV_back_bottom.ini & 
#done

# RIGHT VENTRICLE (back top)
#for SEED in "${SEEDS[@]}"; do
#    ${PROGRAM_PATH} ${INPUT_PATH}/elizabeth_biventricular_coupled_co:length_seed:${SEED}_RV_back_top.ini & 
#done

# RIGHT VENTRICLE (front bottom)
#for SEED in "${SEEDS[@]}"; do
#    ${PROGRAM_PATH} ${INPUT_PATH}/elizabeth_biventricular_coupled_co:length_seed:${SEED}_RV_front_bottom.ini & 
#done

# RIGHT VENTRICLE (front top)
#for SEED in "${SEEDS[@]}"; do
#    ${PROGRAM_PATH} ${INPUT_PATH}/elizabeth_biventricular_coupled_co:length_seed:${SEED}_RV_front_top.ini & 
#done

# LEFT VENTRICLE (copying)
#for SEED in "${SEEDS[@]}"; do
#    cp ${OUTPUT_PATH}/LV/elizabeth_total_length_seed:${SEED}_LV/tree_nterm_520.vtk ${OUTPUT_PATH}/LV/elizabeth_total_length_seed:${SEED}_LV.vtk
#done

# RIGHT VENTRICLE (linking)
#for SEED in "${SEEDS[@]}"; do
#    cp ${OUTPUT_PATH}/RV/elizabeth_total_length_seed:${SEED}_RV_back_bottom/tree_nterm_33.vtk ${TERMINAL_LINKER_PATH}/inputs/elizabeth_min_length_back_bottom.vtk
#    cp ${OUTPUT_PATH}/RV/elizabeth_total_length_seed:${SEED}_RV_back_top/tree_nterm_10.vtk ${TERMINAL_LINKER_PATH}/inputs/elizabeth_min_length_back_top.vtk
#    cp ${OUTPUT_PATH}/RV/elizabeth_total_length_seed:${SEED}_RV_front_bottom/tree_nterm_70.vtk ${TERMINAL_LINKER_PATH}/inputs/elizabeth_min_length_front_bottom.vtk
#    cp ${OUTPUT_PATH}/RV/elizabeth_total_length_seed:${SEED}_RV_front_top/tree_nterm_4.vtk ${TERMINAL_LINKER_PATH}/inputs/elizabeth_min_length_front_top.vtk
#    ${TERMINAL_LINKER_PATH}/bin/TerminalLinker ${TERMINAL_LINKER_PATH}/inputs/elizabeth_min_length_back_top.vtk ${TERMINAL_LINKER_PATH}/inputs/elizabeth_min_length_back_bottom.vtk ${TERMINAL_LINKER_PATH}/inputs/elizabeth_min_length_front_top.vtk ${TERMINAL_LINKER_PATH}/inputs/elizabeth_min_length_front_bottom.vtk ${OUTPUT_PATH}/RV/elizabeth_total_length_seed:${SEED}_RV.vtk
#done

# LVRV (merging)
#for SEED in "${SEEDS[@]}"; do
#    cp ${OUTPUT_PATH}/LV/elizabeth_total_length_seed:${SEED}_LV.vtk ${PURKINJE_MERGER_PATH}/inputs/02_CO_Length_PMJ_last/LV/elizabeth_total_length_seed:${SEED}_LV.vtk
#    cp ${OUTPUT_PATH}/RV/elizabeth_total_length_seed:${SEED}_RV.vtk ${PURKINJE_MERGER_PATH}/inputs/02_CO_Length_PMJ_last/RV/elizabeth_total_length_seed:${SEED}_RV.vtk
#    ${PURKINJE_MERGER_PATH}/bin/Purkinje-Merger ${PURKINJE_MERGER_PATH}/inputs/02_CO_Length_PMJ_last/His/elizabeth_min_length_his.vtk ${PURKINJE_MERGER_PATH}/inputs/02_CO_Length_PMJ_last/LV/elizabeth_total_length_seed:${SEED}_LV.vtk ${PURKINJE_MERGER_PATH}/inputs/02_CO_Length_PMJ_last/RV/elizabeth_total_length_seed:${SEED}_RV.vtk ${PURKINJE_MERGER_PATH}/outputs/02_CO_Length_PMJ_last/elizabeth_total_length_seed:${SEED}_LVRV.vtk
#    cp ${PURKINJE_MERGER_PATH}/outputs/02_CO_Length_PMJ_last/elizabeth_total_length_seed:${SEED}_LVRV.vtk ${MONOALG3D_PATH}/networks/elizabeth_total_length_seed:${SEED}_LVRV.vtk
#done

# ========================================================================================================================
# 3) MINIMIZE TOTAL LENGTH (Multiobjective)
# ========================================================================================================================
#PROGRAM_PATH="/home/berg/Documentos/Elizabeth/Scripts/Network-Generator/bin/Cco_3D"
#INPUT_PATH="/home/berg/Documentos/Elizabeth/Scripts/Network-Generator/inputs/03_CO_Length_Multiobjective"
#OUTPUT_PATH="/home/berg/Documentos/Elizabeth/Scripts/Network-Generator/outputs/03_CO_Length_Multiobjective"
#TERMINAL_LINKER_PATH="/home/berg/Documentos/Elizabeth/Scripts/Terminal-Linker"
#PURKINJE_MERGER_PATH="/home/berg/Documentos/Elizabeth/Scripts/Purkinje-Merger"
#MONOALG3D_PATH="/home/berg/Github/MonoAlg3D_C"

# LEFT VENTRICLE
#for SEED in "${SEEDS[@]}"; do
#    ${PROGRAM_PATH} ${INPUT_PATH}/elizabeth_biventricular_coupled_co:length_seed:${SEED}_LV.ini & 
#done

# RIGHT VENTRICLE (back bottom)
#for SEED in "${SEEDS[@]}"; do
#    ${PROGRAM_PATH} ${INPUT_PATH}/elizabeth_biventricular_coupled_co:length_seed:${SEED}_RV_back_bottom.ini & 
#done

# RIGHT VENTRICLE (back top)
#for SEED in "${SEEDS[@]}"; do
#    ${PROGRAM_PATH} ${INPUT_PATH}/elizabeth_biventricular_coupled_co:length_seed:${SEED}_RV_back_top.ini & 
#done

# RIGHT VENTRICLE (front bottom)
#for SEED in "${SEEDS[@]}"; do
#    ${PROGRAM_PATH} ${INPUT_PATH}/elizabeth_biventricular_coupled_co:length_seed:${SEED}_RV_front_bottom.ini & 
#done

# RIGHT VENTRICLE (front top)
#for SEED in "${SEEDS[@]}"; do
#    ${PROGRAM_PATH} ${INPUT_PATH}/elizabeth_biventricular_coupled_co:length_seed:${SEED}_RV_front_top.ini & 
#done

# LEFT VENTRICLE (copying)
#for SEED in "${SEEDS[@]}"; do
#    cp ${OUTPUT_PATH}/LV/elizabeth_total_length_seed:${SEED}_LV/tree_nterm_520.vtk ${OUTPUT_PATH}/LV/elizabeth_total_length_seed:${SEED}_LV.vtk
#done

# RIGHT VENTRICLE (linking)
#for SEED in "${SEEDS[@]}"; do
#    cp ${OUTPUT_PATH}/RV/elizabeth_total_length_seed:${SEED}_RV_back_bottom/tree_nterm_33.vtk ${TERMINAL_LINKER_PATH}/inputs/elizabeth_min_length_back_bottom.vtk
#    cp ${OUTPUT_PATH}/RV/elizabeth_total_length_seed:${SEED}_RV_back_top/tree_nterm_10.vtk ${TERMINAL_LINKER_PATH}/inputs/elizabeth_min_length_back_top.vtk
#    cp ${OUTPUT_PATH}/RV/elizabeth_total_length_seed:${SEED}_RV_front_bottom/tree_nterm_70.vtk ${TERMINAL_LINKER_PATH}/inputs/elizabeth_min_length_front_bottom.vtk
#    cp ${OUTPUT_PATH}/RV/elizabeth_total_length_seed:${SEED}_RV_front_top/tree_nterm_4.vtk ${TERMINAL_LINKER_PATH}/inputs/elizabeth_min_length_front_top.vtk
#    ${TERMINAL_LINKER_PATH}/bin/TerminalLinker ${TERMINAL_LINKER_PATH}/inputs/elizabeth_min_length_back_top.vtk ${TERMINAL_LINKER_PATH}/inputs/elizabeth_min_length_back_bottom.vtk ${TERMINAL_LINKER_PATH}/inputs/elizabeth_min_length_front_top.vtk ${TERMINAL_LINKER_PATH}/inputs/elizabeth_min_length_front_bottom.vtk ${OUTPUT_PATH}/RV/elizabeth_total_length_seed:${SEED}_RV.vtk
#done

# LVRV (merging)
#for SEED in "${SEEDS[@]}"; do
#    cp ${OUTPUT_PATH}/LV/elizabeth_total_length_seed:${SEED}_LV.vtk ${PURKINJE_MERGER_PATH}/inputs/03_CO_Length_Multiobjective/LV/elizabeth_total_length_seed:${SEED}_LV.vtk
#    cp ${OUTPUT_PATH}/RV/elizabeth_total_length_seed:${SEED}_RV.vtk ${PURKINJE_MERGER_PATH}/inputs/03_CO_Length_Multiobjective/RV/elizabeth_total_length_seed:${SEED}_RV.vtk
#    ${PURKINJE_MERGER_PATH}/bin/Purkinje-Merger ${PURKINJE_MERGER_PATH}/inputs/03_CO_Length_Multiobjective/His/elizabeth_min_length_his.vtk ${PURKINJE_MERGER_PATH}/inputs/03_CO_Length_Multiobjective/LV/elizabeth_total_length_seed:${SEED}_LV.vtk ${PURKINJE_MERGER_PATH}/inputs/03_CO_Length_Multiobjective/RV/elizabeth_total_length_seed:${SEED}_RV.vtk ${PURKINJE_MERGER_PATH}/outputs/03_CO_Length_Multiobjective/elizabeth_total_length_seed:${SEED}_LVRV.vtk
    #cp ${PURKINJE_MERGER_PATH}/outputs/02_CO_Length_PMJ_last/elizabeth_total_length_seed:${SEED}_LVRV.vtk ${MONOALG3D_PATH}/networks/elizabeth_total_length_seed:${SEED}_LVRV.vtk
#done


# ========================================================================================================================
# 4) MINIMIZE TOTAL LENGTH (PMJ's interval)
# ========================================================================================================================
PROGRAM_PATH="/home/berg/Documentos/Elizabeth/Scripts/Network-Generator/bin/Cco_3D"
INPUT_PATH="/home/berg/Documentos/Elizabeth/Scripts/Network-Generator/inputs/04_CO_Length_PMJ_interval_with_region_radius"
OUTPUT_PATH="/home/berg/Documentos/Elizabeth/Scripts/Network-Generator/outputs/04_CO_Length_PMJ_interval_with_region_radius"
TERMINAL_LINKER_PATH="/home/berg/Documentos/Elizabeth/Scripts/Terminal-Linker"
PURKINJE_MERGER_PATH="/home/berg/Documentos/Elizabeth/Scripts/Purkinje-Merger"
MONOALG3D_PATH="/home/berg/Github/MonoAlg3D_C"

# LEFT VENTRICLE
#for SEED in "${SEEDS[@]}"; do
#    ${PROGRAM_PATH} ${INPUT_PATH}/elizabeth_biventricular_coupled_co:length_seed:${SEED}_LV.ini & 
#done

# RIGHT VENTRICLE (back bottom)
#for SEED in "${SEEDS[@]}"; do
#    ${PROGRAM_PATH} ${INPUT_PATH}/elizabeth_biventricular_coupled_co:length_seed:${SEED}_RV_back_bottom.ini & 
#done

# RIGHT VENTRICLE (back top)
#for SEED in "${SEEDS[@]}"; do
#    ${PROGRAM_PATH} ${INPUT_PATH}/elizabeth_biventricular_coupled_co:length_seed:${SEED}_RV_back_top.ini & 
#done

# RIGHT VENTRICLE (front bottom)
#for SEED in "${SEEDS[@]}"; do
#    ${PROGRAM_PATH} ${INPUT_PATH}/elizabeth_biventricular_coupled_co:length_seed:${SEED}_RV_front_bottom.ini & 
#done

# RIGHT VENTRICLE (front top)
#for SEED in "${SEEDS[@]}"; do
#    ${PROGRAM_PATH} ${INPUT_PATH}/elizabeth_biventricular_coupled_co:length_seed:${SEED}_RV_front_top.ini & 
#done

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
for SEED in "${SEEDS[@]}"; do
    cp ${OUTPUT_PATH}/LV/elizabeth_total_length_seed:${SEED}_LV.vtk ${PURKINJE_MERGER_PATH}/inputs/04_CO_Length_PMJ_interval_with_region_radius/LV/elizabeth_total_length_seed:${SEED}_LV.vtk
    cp ${OUTPUT_PATH}/RV/elizabeth_total_length_seed:${SEED}_RV.vtk ${PURKINJE_MERGER_PATH}/inputs/04_CO_Length_PMJ_interval_with_region_radius/RV/elizabeth_total_length_seed:${SEED}_RV.vtk
    ${PURKINJE_MERGER_PATH}/bin/Purkinje-Merger ${PURKINJE_MERGER_PATH}/inputs/04_CO_Length_PMJ_interval_with_region_radius/His/elizabeth_min_length_his.vtk ${PURKINJE_MERGER_PATH}/inputs/04_CO_Length_PMJ_interval_with_region_radius/LV/elizabeth_total_length_seed:${SEED}_LV.vtk ${PURKINJE_MERGER_PATH}/inputs/04_CO_Length_PMJ_interval_with_region_radius/RV/elizabeth_total_length_seed:${SEED}_RV.vtk ${PURKINJE_MERGER_PATH}/outputs/04_CO_Length_PMJ_interval_with_region_radius/elizabeth_total_length_seed:${SEED}_LVRV.vtk
    #cp ${PURKINJE_MERGER_PATH}/outputs/02_CO_Length_PMJ_last/elizabeth_total_length_seed:${SEED}_LVRV.vtk ${MONOALG3D_PATH}/networks/elizabeth_total_length_seed:${SEED}_LVRV.vtk
done