#!/bin/bash

run_simulation_Noble () {
    echo "****************************************************************************************************"
    echo "[!] Running Noble simulation with $1ms pacing ..."
    echo "****************************************************************************************************"
    # Steady-State 
    ./bin/FDMMonodomain1D examples/sst_cable_5cm_$1ms.ini
    cp output.sst steady_state/cable-5cm-$1ms.sst 
    #./clear_results.sh
    # Experiment
    ./bin/FDMMonodomain1D examples/simple_cable_5cm_$1ms.ini
    mkdir scripts/$1ms
    cp output/*.dat scripts/$1ms
    echo "****************************************************************************************************"
}


#PACINGS=( 200 210 220 230 240 250 260 270 280 290 300 )
PACINGS=($(seq 280 -5 200))

for PACING in "${PACINGS[@]}"; do
    run_simulation_Noble $PACING
done
