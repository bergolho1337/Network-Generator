seeds = [1562002891,1562002894,1562005513,1562005553,1562006177,1562007596,1562008172,1562008424,1562009134,1562009769]

# RV
for seed in seeds:
    file = open("elizabeth_biventricular_coupled_co:length_seed:%d_RV.ini" % (seed),"w")
    file.write('''[main]
N_term = 496
root_x = 0.05275
root_y = 0.0605	
root_z = 0.0185
max_rand_offset = 2
seed = %d
use_only_murray = false

[save_network]
output_dir = outputs/04_Elizabeth/01_CO_Length/RV_seed:%d

[cloud_points]
use_cloud_points = true
cloud_points_filename = clouds/private/01_Elizabeth_Cherry/elizabeth_cloud_remapped_496term.pts
use_obstacle = false
use_pmj_location = false

[local_optimization]
use_local_optimization = true
local_optimization_function = rafael_local_optimization

[cost_function]
library_name = shared-libs/libminimize_custom_function.so
function_name = minimize_custom_function_with_angle_restriction
beta = 1.0
alpha = 0.0
min_degrees_limit = 1
max_degrees_limit = 85
''' % (seed,seed))
    file.close()

# LV
for seed in seeds:
    file = open("elizabeth_biventricular_coupled_co:length_seed:%d_LV.ini" % (seed),"w")
    file.write('''[main]
N_term = 650
root_x = 0.04  
root_y = 0.061
root_z = 0.010
max_rand_offset = 2
seed = %d
use_only_murray = false

[save_network]
output_dir = outputs/04_Elizabeth/01_CO_Length/LV_seed:%d

[cloud_points]
use_cloud_points = true
cloud_points_filename = clouds/private/01_Elizabeth_Cherry/elizabeth_cloud_remapped_650term_LV.pts
use_obstacle = false
use_pmj_location = false

[local_optimization]
use_local_optimization = true
local_optimization_function = rafael_local_optimization

[cost_function]
library_name = shared-libs/libminimize_custom_function.so
function_name = minimize_custom_function_with_angle_restriction
beta = 1.0
alpha = 0.0
min_degrees_limit = 1
max_degrees_limit = 63
''' % (seed,seed))
    file.close()