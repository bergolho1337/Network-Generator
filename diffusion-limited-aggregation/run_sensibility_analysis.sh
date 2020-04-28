# Author: Lucas Berg
# ------------------------------------------------------------------------------------------------------------ 
# This script is able to generate a set Purkinje networks, which the configuration files should be previously 
# configured with the sensibility-analysis script.
# ------------------------------------------------------------------------------------------------------------ 
#!/bin/bash

PROGRAM_NAME="./bin/DLA"
INPUT_FOLDER="inputs/dla_config"
TRANSFORM_SCRIPT_NAME="./scripts/transform-filter/bin/TransformFilter"
MONOALG_NETWORK_FOLDER="/home/berg/Github/MonoAlg3D_C/networks/elizabeth-meshes/dla"

SEEDS=( 1562046115 1562013988 1562042299 1562005512 1562009134 1562009768 1562044566 1562008423 1562036996 1562020974 )
WALKER_RADIUS=( 2.1 2.2 2.3 2.4 2.5 2.6 2.7 2.8 2.9 )
COUNTER=1

for RADIUS in "${WALKER_RADIUS[@]}"; do 
    for SEED in "${SEEDS[@]}"; do
            
        # Run the Purkinje network generator
        ${PROGRAM_NAME} "${INPUT_FOLDER}/elizabeth_purkinje_dla_seed-${SEED}_radius-${RADIUS}.ini"
        
        # Move the output to the pos-processing folder
        mv output/dla_tree.vtk output/dla_chaos/elizabeth_purkinje_dla_${COUNTER}.vtk 

        # Pos-process the network by parsing from the CCO domain to the MonoAlg3D one.
        ${TRANSFORM_SCRIPT_NAME} output/dla_chaos/elizabeth_purkinje_dla_${COUNTER}.vtk ${MONOALG_NETWORK_FOLDER}/elizabeth_purkinje_dla_seed-${SEED}_radius-${RADIUS}.vtk

        let "COUNTER=COUNTER+1"
    done
done

# After this script the Purkinje network will be already into the MonoAlg3D domain and ready to be activated by the Monodomain solver.
# Just move the monoalg3d configuration files to the MonoAlg3D_C/example_configs/ folder and run the batch script.

