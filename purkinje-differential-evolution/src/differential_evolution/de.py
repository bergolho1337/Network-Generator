import os
import subprocess
import network_helper as net
import monoalg_helper as mono
import numpy as np
from scipy.optimize import differential_evolution

SEEDS=[1562002891,1562002894,1562005513,1562005553,1562006177,1562007596,1562008172,1562008424,1562009135,1562009769]
#SEEDS=[1562002891,1562002894,1562005513,1562005553,1562006177]

# =============================================================================================================================
# Objective function for the DE
def purkinje_function (x):
    # We can have multiple parameters to fit here 
    C = x[0]        # Pruning function parameter
    
    net.generate_network_generator_config_files(SEEDS,C)
    net.run_network_generator(SEEDS,C)
    
    mono.generate_monoalg3d_config_files(SEEDS)
    rrmse_array = mono.run_monoalg3d(SEEDS)

    result = np.mean(rrmse_array)
    print(rrmse_array)
    print("\tC = %g || Eval = %g" % (C,result))
    return result

# =============================================================================================================================

def create_folders ():
    for seed in SEEDS:
        subprocess.call(["mkdir","../outputs/networks/seed:%d" % (seed)])
        subprocess.call(["mkdir","../outputs/lats/seed:%d" % (seed)])
        subprocess.call(["mkdir","../outputs/errors/seed:%d" % (seed)])

def solve_differential_evolution ():

    print("[solve_differential_evolution] Creating folders ...")
    create_folders()

    # Define the objective function and the bounds
    fn = purkinje_function
    bounds = [(0,5)]

    result = differential_evolution(fn,bounds,maxiter=5,popsize=5,disp=True)

    if (result.success):
        print("\tSuccess!")
    else:
        print("\tERROR!")
    print("\tSolution = [%s]" % (result.x))
    print("\tMessage = %s" % (result.message))
    print("\tFunction value = %s" % (result.fun))
    print("\tNumber of objective function evaluations = %s" % (result.nfev))
    print("\tNumber of iterations = %s" % (result.nit))