# Author: Lucas Berg
# ------------------------------------------------------------------------------------------------------------ 
# This script is able to generate a set Purkinje networks. The configuration files should be previously 
# configured with the sensibility-analysis script.
# ------------------------------------------------------------------------------------------------------------ 
#!/bin/bash

# HELPER FUNCTIONS
function run_co_algorithm_tests () {

	PROGRAM_NAME=${1}
	TRANSFORM_SCRIPT_NAME=${2}
	SEEDS=${3}
	RAND_OFFSETS=${4}

	INPUT_FOLDER="inputs/co_length_config"
	MONOALG_NETWORK_FOLDER="/home/berg/Github/MonoAlg3D_C/networks/elizabeth-meshes/co_length"

	COUNTER=1
	for RAND_OFFSET in "${RAND_OFFSETS[@]}"; do
    		for SEED in "${SEEDS[@]}"; do    
            		#echo "elizabeth_purkinje_cco_seed-${SEED}_offset-${RAND_OFFSET}_nterm-130.ini"

        		# Run the Purkinje network generator
        		${PROGRAM_NAME} "${INPUT_FOLDER}/elizabeth_purkinje_co_seed-${SEED}_offset-${RAND_OFFSET}_nterm-130.ini"

        		# Move the output to the pos-processing folder
        		mv output/cco_tree_cm.vtk output/elizabeth_chaos_co/elizabeth_purkinje_co_${COUNTER}.vtk 

        		# Pos-process the network by parsing from the CCO domain to the MonoAlg3D one.
        		${TRANSFORM_SCRIPT_NAME} output/elizabeth_chaos_co/elizabeth_purkinje_co_${COUNTER}.vtk ${MONOALG_NETWORK_FOLDER}/elizabeth_purkinje_co_seed-${SEED}_offset-${RAND_OFFSET}_nterm-130.vtk

        		let "COUNTER=COUNTER+1"
    		done
	done

}

function run_co_activation_time_algorithm_tests () {

	PROGRAM_NAME=${1}
	TRANSFORM_SCRIPT_NAME=${2}
	SEEDS=${3}
	RAND_OFFSETS=${4}

	INPUT_FOLDER="inputs/co_activation_time_config"
	MONOALG_NETWORK_FOLDER="/home/berg/Github/MonoAlg3D_C/networks/elizabeth-meshes/co_activation_time"

	COUNTER=1
	for RAND_OFFSET in "${RAND_OFFSETS[@]}"; do
    		for SEED in "${SEEDS[@]}"; do    
            		#echo "elizabeth_purkinje_cco_seed-${SEED}_offset-${RAND_OFFSET}_nterm-130.ini"

        		# Run the Purkinje network generator
        		${PROGRAM_NAME} "${INPUT_FOLDER}/elizabeth_purkinje_co_activation_time_seed-${SEED}_offset-${RAND_OFFSET}_nterm-130.ini"

        		# Move the output to the pos-processing folder
        		mv output/cco_tree_cm.vtk output/elizabeth_chaos_co_activation_time/elizabeth_purkinje_co_${COUNTER}.vtk 

        		# Pos-process the network by parsing from the CCO domain to the MonoAlg3D one.
        		${TRANSFORM_SCRIPT_NAME} output/elizabeth_chaos_co_activation_time/elizabeth_purkinje_co_${COUNTER}.vtk ${MONOALG_NETWORK_FOLDER}/elizabeth_purkinje_co_activation_time_seed-${SEED}_offset-${RAND_OFFSET}_nterm-130.vtk

        		let "COUNTER=COUNTER+1"
    		done
	done

}

function run_cco_algorithm_tests () {

        PROGRAM_NAME=${1}
        TRANSFORM_SCRIPT_NAME=${2}
        SEEDS=${3}
        RAND_OFFSETS=${4}

		INPUT_FOLDER="inputs/cco_classic_config"
        MONOALG_NETWORK_FOLDER="/home/berg/Github/MonoAlg3D_C/networks/elizabeth-meshes/cco_classic"

        COUNTER=1
        for RAND_OFFSET in "${RAND_OFFSETS[@]}"; do
                for SEED in "${SEEDS[@]}"; do    
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

}


# MAIN PROGRAM
PROGRAM_NAME="./bin/Cco_3D"
TRANSFORM_SCRIPT_NAME="./scripts/transform-filter/bin/TransformFilter"

SEEDS=( 1562046115 1562013988 1562042299 1562005513 1562009134 1562009769 1562044567 1562008424 1562036996 1562020974 1562049673 1562023024 1562025823 1562007596 1562028813 1562005553 1562017900 1562034865 1562047507 1562011558 1562051289 1562006177 1562002894 1562002891 1562018879 1562021648 1562035947 1562008172 1562043344 1562021993 )

#RAND_OFFSETS=( 2 3 4 5 6 7 8 9 10 )
RAND_OFFSETS=( 9 )


# RUNNING ALL THE ALGORITHM TESTS ...
#run_co_algorithm_tests $PROGRAM_NAME $TRANSFORM_SCRIPT_NAME $SEEDS $RAND_OFFSETS
run_co_activation_time_algorithm_tests $PROGRAM_NAME $TRANSFORM_SCRIPT_NAME $SEEDS $RAND_OFFSETS
#run_cco_algorithm_tests $PROGRAM_NAME $TRANSFORM_SCRIPT_NAME $SEEDS $RAND_OFFSETS

# After this script the Purkinje network will be already into the MonoAlg3D domain and ready to be activated by the Monodomain solver.
# Just move the monoalg3d configuration files to the MonoAlg3D_C/example_configs/ folder and run the batch script.

