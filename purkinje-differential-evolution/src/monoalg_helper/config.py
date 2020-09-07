import os
import subprocess
import numpy as np

MONOALG3D_CONFIG_FILE_PATH = "/home/berg/Github/Network-Generator/purkinje-differential-evolution/monoalg-config"
MONOALG3D_OUTPUT_DIR_PATH = "/home/berg/Github/Network-Generator/purkinje-differential-evolution/outputs/lats"
MONOALG3D_SHARED_LIBS = "/home/berg/Github/MonoAlg3D_C/shared_libs"
MONOALG3D_BINARY_PATH = "/home/berg/Github/MonoAlg3D_C/bin"
MONOALG3D_PMJ_LOCATION_PATH = "/home/berg/Github/MonoAlg3D_C/networks/03_Lucas/02_Oxford/01_Basics/00_Reference"

ERROR_CALCULATOR_REFERENCE_LAT_PATH = "/home/berg/Github/MonoAlg3D_C/scripts/00_Tools/01_Error_Calculator/scripts/03_Oxford/inputs/00_Reference/activation_map_gold.vtu"
ERROR_CALCULATOR_OUTPUT_PATH = "/home/berg/Github/Network-Generator/purkinje-differential-evolution/outputs/errors"
ERROR_CALCULATOR_BINARY_PATH = "/home/berg/Github/MonoAlg3D_C/scripts/00_Tools/01_Error_Calculator/bin"

NETWORK_GENERATOR_OUTPUT_DIR_PATH = "/home/berg/Github/Network-Generator/purkinje-differential-evolution/outputs/networks"

DEVNULL = open(os.devnull, 'w')

def generate_monoalg3d_config_files (seeds):
    for seed in seeds:
        write_monoalg3d_config_file(seed)

def write_monoalg3d_config_file (seed):
    filename = "%s/seed_%d.ini" % (MONOALG3D_CONFIG_FILE_PATH,seed)

    file = open(filename,"w")
    file.write('''[main]
num_threads = 6
dt_pde = 0.02
simulation_time = 150.0
abort_on_no_activity = false
use_adaptivity = false

[update_monodomain]
main_function = update_monodomain_default
library_file = %s/libdefault_update_monodomain.so

[save_result]
print_rate = 10
output_dir = %s/seed:%d
init_function=init_save_purkinje_coupling_with_activation_times
end_function=end_save_purkinje_coupling_with_activation_times
main_function=save_purkinje_coupling_with_activation_times
activation_threshold_tissue = 0.3
activation_threshold_purkinje = 0.3
apd_threshold_tissue = 0.1
apd_threshold_purkinje = 0.1
save_activation_time = true
save_apd = false
save_pvd = true
file_prefix = V_Tissue
file_prefix_purkinje = V_Purkinje
remove_older_simulation=true
library_file = %s/libdefault_save_mesh.so

[assembly_matrix]
init_function = set_initial_conditions_coupling_fvm
sigma_x = 0.0000176
sigma_y = 0.0000176
sigma_z = 0.0000176
sigma_purkinje = 0.001
main_function = purkinje_coupling_assembly_matrix
library_file = %s/libpurkinje_coupling_matrix_assembly.so

[linear_system_solver]
tolerance = 1e-16
use_preconditioner = no
max_iterations = 200
main_function = conjugate_gradient
library_file = %s/libdefault_linear_system_solver.so

[alg]
refinement_bound = 0.11
derefinement_bound = 0.10
refine_each = 1
derefine_each = 1

[domain]
name = Plain Mesh
num_layers = 1
start_dx = 100.0
start_dy = 100.0
start_dz = 100.0
side_length=20000
main_function = initialize_grid_with_square_mesh
library_file = %s/libdefault_domains.so

[purkinje]
name = Purkinje
start_discretization = 100.0
rpmj = 1000.0
pmj_scale = 500.0
nmin_pmj = 10
nmax_pmj = 30
retro_propagation = true
network_file = %s/seed:%d/converted_network.vtk
pmj_location_file = %s/reference_pmjs.vtk
main_function = initialize_purkinje_coupling_with_custom_mesh
library_file = %s/libdefault_purkinje.so

[purkinje_ode_solver]
dt = 0.02
use_gpu = no
gpu_id = 0
library_file = %s/libfhn_mod.so

[ode_solver]
dt = 0.02
use_gpu = yes
gpu_id = 0
library_file = %s/libmitchell_shaeffer_2003.so

[purkinje_stim_his]
start = 0.0
duration = 2.0
current = 1.0
x_limit = 500
main_function = stim_if_x_less_than
library_file = %s/libdefault_stimuli.so

               ''' % (MONOALG3D_SHARED_LIBS,MONOALG3D_OUTPUT_DIR_PATH,seed,MONOALG3D_SHARED_LIBS,MONOALG3D_SHARED_LIBS,MONOALG3D_SHARED_LIBS,MONOALG3D_SHARED_LIBS,NETWORK_GENERATOR_OUTPUT_DIR_PATH,seed,MONOALG3D_PMJ_LOCATION_PATH,MONOALG3D_SHARED_LIBS,MONOALG3D_SHARED_LIBS,MONOALG3D_SHARED_LIBS,MONOALG3D_SHARED_LIBS))
    file.close()

def run_monoalg3d (seeds):
    # For each seed run the MonoAlg3D solver and calculate its RRMSE error
    for seed in seeds:
        # Call the MonoAlg3D solver supressing the stdout output
        subprocess.call(["%s/MonoAlg3D" % (MONOALG3D_BINARY_PATH),"-c","%s/seed_%d.ini" % (MONOALG3D_CONFIG_FILE_PATH,seed)], stdout=DEVNULL)

        # Call the ErrorCalculator app to convert the Purkinje network to the MonoAlg3D domain
        STDOUT = open("%s/seed:%d/rrmse.dat" % (ERROR_CALCULATOR_OUTPUT_PATH,seed),"w")
        subprocess.call(["%s/ErrorCalculator" % (ERROR_CALCULATOR_BINARY_PATH),"%s" % (ERROR_CALCULATOR_REFERENCE_LAT_PATH),"%s/seed:%d/tissue_activation_time_map_pulse_it_0.vtu" % (MONOALG3D_OUTPUT_DIR_PATH,seed),"%s/seed:%d/error.vtu" % (ERROR_CALCULATOR_OUTPUT_PATH,seed)], stdout=STDOUT)

    # Now open each RRMSE file and build an array with the results from the Purkinje networks
    rrmse_array = []
    for seed in seeds:
        rrmse = np.genfromtxt("%s/seed:%d/rrmse.dat" % (ERROR_CALCULATOR_OUTPUT_PATH,seed))
        rrmse_array.append(rrmse)
    
    return rrmse_array
