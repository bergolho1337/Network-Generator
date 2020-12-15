import sys
import numpy as np

def print_config_LV (seeds):

    for seed in seeds:
        file = open("elizabeth_biventricular_coupled_co:length_seed:%d_LV.ini" % (seed),"w")
        file.write('''[main]
N_term = 520
root_x = 0.0344435
root_y = 0.0554435
root_z = 0.00644345
max_rand_offset = 2
seed = %d
use_only_murray = true
start_radius = 0.1
gamma = 19

[save_network]
output_dir = outputs/04_PMJ_Linktries:40/LV_seed:%d

[cloud_points]
use_cloud_points = true
cloud_points_filename = clouds/inactives/jesuliana_mesh_cloud_LV_m.vtk
use_pmj_location = true
pmj_location_filename = clouds/pmjs/jesuliana_mesh_active_pmjs_LV_m.vtk
max_pmj_connection_tries = 50
pmj_connection_rate = 13
pmj_region_radius = 0.002
lat_error_tolerance = 2.0
lat_offset = 5.2135
use_obstacle = false

[local_optimization]
use_local_optimization = true
local_optimization_function = rafael_local_optimization

[cost_function]
library_name = shared-libs/libminimize_custom_function.so
function_name = minimize_custom_function
beta = 1.0
alpha = 0.0
min_degrees_limit = 1
max_degrees_limit = 63
min_segment_length = 0.0001
max_segment_length = 0.007
    ''' % (seed,seed))
        file.close()

def print_config_back_top (seeds):
    for seed in seeds:
        file = open("elizabeth_biventricular_coupled_co:length_seed:%d_RV_back_top.ini" % (seed),"w")
        file.write('''[main]
N_term = 10
root_x = 0.0469435
root_y = 0.0571935
root_z = 0.0126935
max_rand_offset = 2
seed = %d
use_only_murray = true
start_radius = 0.1
gamma = 19

[save_network]
output_dir = outputs/04_PMJ_Linktries:40/RV_back_top_seed:%d

[cloud_points]
use_cloud_points = true
cloud_points_filename = clouds/inactives/jesuliana_mesh_cloud_RV_back_top_m.vtk
use_pmj_location = false
use_obstacle = false
lat_offset = 7.70042

[local_optimization]
use_local_optimization = true
local_optimization_function = rafael_local_optimization

[cost_function]
library_name = shared-libs/libminimize_custom_function.so
function_name = minimize_custom_function
beta = 1.0
alpha = 0.0
min_degrees_limit = 1
max_degrees_limit = 65
min_segment_length = 0.0001
max_segment_length = 0.01    
    ''' % (seed,seed))
        file.close()

def print_config_back_bottom (seeds):
    for seed in seeds:
        file = open("elizabeth_biventricular_coupled_co:length_seed:%d_RV_back_bottom.ini" % (seed),"w")
        file.write('''[main]
N_term = 40
root_x = 0.0579435
root_y = 0.0506935
root_z = 0.0258835
max_rand_offset = 2
seed = %d
use_only_murray = true
start_radius = 0.1
gamma = 19

[save_network]
output_dir = outputs/04_PMJ_Linktries:40/RV_back_bottom_seed:%d

[cloud_points]
use_cloud_points = true
cloud_points_filename = clouds/inactives/jesuliana_mesh_cloud_RV_back_bottom_m.vtk
use_pmj_location = true
pmj_location_filename = clouds/pmjs/jesuliana_mesh_active_pmjs_RV_back_bottom_m.vtk
max_pmj_connection_tries = 200
pmj_connection_rate = 1
pmj_region_radius = 0.002
lat_error_tolerance = 2.0
use_obstacle = false

[local_optimization]
use_local_optimization = true
local_optimization_function = rafael_local_optimization

[cost_function]
library_name = shared-libs/libminimize_custom_function.so
function_name = minimize_custom_function
beta = 1.0
alpha = 0.0
min_degrees_limit = 1
max_degrees_limit = 65
min_segment_length = 0.0001
max_segment_length = 0.01
    ''' % (seed,seed))
        file.close()

def print_config_front_top (seeds):
    for seed in seeds:
        file = open("elizabeth_biventricular_coupled_co:length_seed:%d_RV_front_top.ini" % (seed),"w")
        file.write('''[main]
N_term = 4
root_x = 0.0559435
root_y = 0.0489552
root_z = 0.0126935
max_rand_offset = 2
seed = %d
use_only_murray = true
start_radius = 0.1
gamma = 19

[save_network]
output_dir = outputs/04_PMJ_Linktries:40/RV_front_top_seed:%d

[cloud_points]
use_cloud_points = true
cloud_points_filename = clouds/inactives/jesuliana_mesh_cloud_RV_front_top_m.vtk
use_pmj_location = true
pmj_location_filename = clouds/pmjs/jesuliana_mesh_active_pmjs_RV_front_top_m.vtk
max_pmj_connection_tries = 200 
pmj_connection_rate = 4
pmj_region_radius = 0.002
lat_error_tolerance = 2.0
use_obstacle = false

[local_optimization]
use_local_optimization = true
local_optimization_function = rafael_local_optimization

[cost_function]
library_name = shared-libs/libminimize_custom_function.so
function_name = minimize_custom_function
beta = 1.0
alpha = 0.0
min_degrees_limit = 1
max_degrees_limit = 65
min_segment_length = 0.0001
max_segment_length = 0.01
    ''' % (seed,seed))
        file.close()

def print_config_front_bottom (seeds):
    for seed in seeds:
        file = open("elizabeth_biventricular_coupled_co:length_seed:%d_RV_front_bottom.ini" % (seed),"w")
        file.write('''[main]
N_term = 40
root_x = 0.0577575
root_y = 0.0484435
root_z = 0.0171934
max_rand_offset = 2
seed = %d
use_only_murray = true
start_radius = 0.1
gamma = 19

[save_network]
output_dir = outputs/04_PMJ_Linktries:40/RV_front_bottom_seed:%d

[cloud_points]
use_cloud_points = true
cloud_points_filename = clouds/inactives/jesuliana_mesh_cloud_RV_front_bottom_m.vtk
use_pmj_location = true
pmj_location_filename = clouds/pmjs/jesuliana_mesh_active_pmjs_RV_front_bottom_m.vtk
max_pmj_connection_tries = 200
pmj_connection_rate = 1
pmj_region_radius = 0.002
lat_error_tolerance = 2.0
use_obstacle = false

[local_optimization]
use_local_optimization = true
local_optimization_function = rafael_local_optimization

[cost_function]
library_name = shared-libs/libminimize_custom_function.so
function_name = minimize_custom_function
beta = 1.0
alpha = 0.0
min_degrees_limit = 1
max_degrees_limit = 65
min_segment_length = 0.0001
max_segment_length = 0.01
    ''' % (seed,seed))
        file.close()

def print_config_RV (seeds):
    print_config_back_top(seeds)
    print_config_back_bottom(seeds)
    print_config_front_top(seeds)
    print_config_front_bottom(seeds)


def main():
    
    #seeds = [1562002891,1562002894,1562005513,1562005553,1562006177,1562007596,1562008172,1562008424,1562009134,1562009769]
    seeds = [1562002891,1562002894,1562005513,1562005555,1562006177]
    
    print_config_LV(seeds)
    print_config_RV(seeds)

if __name__ == "__main__":
        main()
