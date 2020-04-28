# Author: Lucas Berg
# ------------------------------------------------------------------------------------------------------------ 
# This script is able to generate a set Purkinje networks, which the configuration files should be previously 
# configured with the sensibility-analysis script.
# ------------------------------------------------------------------------------------------------------------ 
#!/bin/bash

PROGRAM_NAME="./bin/Cco_3D"
INPUT_FOLDER="inputs/cco_config"
TRANSFORM_SCRIPT_NAME="./scripts/transform-filter/bin/TransformFilter"
MONOALG_NETWORK_FOLDER="/home/berg/Github/MonoAlg3D_C/networks/elizabeth-meshes/cco_classic"

SEEDS=( 1562046115 1562013988 1562042299 1562005512 1562009134 1562009768 1562044566 1562008423 1562036996 1562020974 )
#SEEDS=( 1562046115 )
RAND_OFFSETS=( 2 3 4 5 6 7 8 9 10 )
#RAND_OFFSETS=( 6 )
COUNTER=1

for SEED in "${SEEDS[@]}"; do
    for RAND_OFFSET in "${RAND_OFFSETS[@]}"; do
	#echo "elizabeth_purkinje_cco_seed-${SEED}_offset-${RAND_OFFSET}_nterm-130.ini"
        
        # Run the Purkinje network generator
        ${PROGRAM_NAME} "${INPUT_FOLDER}/elizabeth_purkinje_cco_seed-${SEED}_offset-${RAND_OFFSET}_nterm-130.ini"
        
        # Move the output to the pos-processing folder
        mv output/cco_tree_cm.vtk output/elizabeth_chaos_cco/elizabeth_purkinje_cco_${COUNTER}.vtk 

        # Pos-process the network by parsing from the CCO domain to the MonoAlg3D one.
        ${TRANSFORM_SCRIPT_NAME} output/elizabeth_chaos_cco/elizabeth_purkinje_cco_${COUNTER}.vtk ${MONOALG_NETWORK_FOLDER}/elizabeth_purkinje_cco_seed-${SEED}_offset-${RAND_OFFSET}_nterm-130.vtk

        let "COUNTER=COUNTER+1"
    done
done

# After this script the Purkinje network will be already into the MonoAlg3D domain and ready to be activated by the Monodomain solver.
# Just move the monoalg3d configuration files to the MonoAlg3D_C/example_configs/ folder and run the batch script.

