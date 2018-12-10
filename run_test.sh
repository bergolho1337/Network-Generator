#!/bin/bash

run_test_Noble () {
    echo "****************************************************************************************************"
    echo "[!] Running Noble test ..."
    echo "****************************************************************************************************"
    # Steady-State 
    ./bin/FDMMonodomain1D examples/sst_sample.ini
    cp output.sst steady_state/sst_sample.sst 
    #./clear_results.sh
    # Experiment
    ./bin/FDMMonodomain1D examples/simple_sample.ini
    echo "****************************************************************************************************"
}

run_test_Noble

