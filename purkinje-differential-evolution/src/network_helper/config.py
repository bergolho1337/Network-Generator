import os
import subprocess

NETWORK_GENERATOR_CONFIG_FILE_PATH = "/home/berg/Github/Network-Generator/purkinje-differential-evolution/network-config"
NETWORK_GENERATOR_OUTPUT_DIR_PATH = "/home/berg/Github/Network-Generator/purkinje-differential-evolution/outputs/networks"
NETWORK_GENERATOR_SHARED_LIBS = "/home/berg/Github/Network-Generator/cco-project-3d/shared-libs"
NETWORK_GENERATOR_BINARY_PATH = "/home/berg/Github/Network-Generator/cco-project-3d/bin"
NETWORK_GENERATOR_CLOUD_POINTS_PATH = "/home/berg/Github/Network-Generator/cco-project-3d/clouds/generated-clouds"

TRANSFORM_FILTER_BINARY_PATH = "/home/berg/Github/Network-Generator/cco-project-3d/scripts/transform-filter/bin"

DEVNULL = open(os.devnull, 'w')

def generate_network_generator_config_files (seeds, C):
    for seed in seeds:
        write_network_generator_config_file(seed,C)

def write_network_generator_config_file (seed, C):
    filename = "%s/seed_%d.ini" % (NETWORK_GENERATOR_CONFIG_FILE_PATH,seed)

    file = open(filename,"w")
    file.write('''[main]
Q_perf = 8.33e-06
p_perf = 13300
p_term = 9580
V_perf = 0.0001
N_term = 80
root_x = 0.001
root_y = 0.0985
root_z = 0.0005
max_rand_offset = 2
seed = %d
use_only_murray = true
start_radius = 0.1
gamma = 19

[save_network]
output_dir = %s/seed:%d

[cloud_points]
use_cloud_points = true
cloud_points_filename = %s/slab_coordinates_guided.pts
use_obstacle = false
use_pmj_location = true
pmj_location_filename = %s/slab_pmj_locations.pts

[local_optimization]
use_local_optimization = true
local_optimization_function = rafael_local_optimization
library_name = %s/libdefault_local_optimization.so

[cost_function]
library_name = %s/libminimize_activation_time.so
function_name = minimize_weighted_activation_time_with_angle_restriction
w = 1
G = 7.9
Cf = 3.4
tauf = 0.1
min_degrees_limit = 1
max_degrees_limit = 63

[pruning]
use_pruning = true
library_name = %s/libdefault_pruning.so
pruning_function = hyperbolic_tangent
A = 50.0
B = -0.25
C = %g
D = 50.0
               ''' % (seed,NETWORK_GENERATOR_OUTPUT_DIR_PATH,seed,NETWORK_GENERATOR_CLOUD_POINTS_PATH,NETWORK_GENERATOR_CLOUD_POINTS_PATH,NETWORK_GENERATOR_SHARED_LIBS,NETWORK_GENERATOR_SHARED_LIBS,NETWORK_GENERATOR_SHARED_LIBS,C))
    file.close()

def run_network_generator (seeds,C):
    # For each seed generate the Purkinje network using the configuration file and them apply a transformation filter to convert to the MonoAlg3D domain
    for seed in seeds:
        # Call the Network-Generator solver supressing the stdout output
        subprocess.call(["%s/Cco_3D" % (NETWORK_GENERATOR_BINARY_PATH),"%s/seed_%d.ini" % (NETWORK_GENERATOR_CONFIG_FILE_PATH,seed)], stdout=DEVNULL)

        # Call the Transform-Filter app to convert the Purkinje network to the MonoAlg3D domain
        subprocess.call(["%s/TransformFilter" % (TRANSFORM_FILTER_BINARY_PATH),"%s/seed:%d/tree_nterm_80.vtk" % (NETWORK_GENERATOR_OUTPUT_DIR_PATH,seed), "%s/seed:%d/converted_network.vtk" % (NETWORK_GENERATOR_OUTPUT_DIR_PATH,seed)])
