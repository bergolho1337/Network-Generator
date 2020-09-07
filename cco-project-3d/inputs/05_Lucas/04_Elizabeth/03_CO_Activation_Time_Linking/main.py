seeds = [1562002891,1562002894,1562005513,1562005553,1562006177,1562007596,1562008172,1562008424,1562009134,1562009769]

# RV
for seed in seeds:
    file = open("elizabeth_biventricular_coupled_co:at_seed:%d_RV.ini" % (seed),"w")
    file.write('''[main]
Q_perf = 8.33e-06
p_perf = 13300
p_term = 9580
p_perf = 13300
V_perf = 0.0001
N_term = 140
root_x = 0.051
root_y = 0.062453
root_z = 0.01925
max_rand_offset = 2
seed = %d
use_only_murray = true
start_radius = 0.1
gamma = 19

[save_network]
output_dir = outputs/04_Elizabeth/03_CO_Activation_Time_Linking/RV_seed:%d

[cloud_points]
use_cloud_points = true
cloud_points_filename = clouds/private/01_Elizabeth_Cherry/elizabeth_guided_cloud_RV_r:0.007.pts
use_obstacle = false
use_pmj_location = true
pmj_location_filename = clouds/private/01_Elizabeth_Cherry/elizabeth_pmjs_rescale_RV.pts

[local_optimization]
use_local_optimization = true
local_optimization_function = rafael_local_optimization

[cost_function]
library_name = shared-libs/libminimize_activation_time.so
function_name = minimize_activation_time_with_angle_restriction
G = 7.9
Cf = 3.4
tauf = 0.1
min_degrees_limit = 1
max_degrees_limit = 70
''' % (seed,seed))
    file.close()

# LV
for seed in seeds:
    file = open("elizabeth_biventricular_coupled_co:at_seed:%d_LV.ini" % (seed),"w")
    file.write('''[main]
Q_perf = 8.33e-06
p_perf = 13300
p_term = 9580
p_perf = 13300
V_perf = 0.0001
N_term = 650
root_x = -0.004
root_y = 0.02
root_z = 0.0265
max_rand_offset = 2
seed = %d
use_only_murray = true
start_radius = 0.1
gamma = 19

[save_network]
output_dir = outputs/04_Elizabeth/03_CO_Activation_Time_Linking/LV_seed:%d

[cloud_points]
use_cloud_points = true
cloud_points_filename = clouds/private/elizabeth_remapped_guided_3.pts
use_obstacle = false
use_pmj_location = true
pmj_location_filename = clouds/private/elizabeth_exterior_LV_pmjs.pts

[local_optimization]
use_local_optimization = true
local_optimization_function = rafael_local_optimization

[cost_function]
library_name = shared-libs/libminimize_activation_time.so
function_name = minimize_activation_time_with_angle_restriction
G = 7.9
Cf = 3.4
tauf = 0.1
min_degrees_limit = 1
max_degrees_limit = 63
''' % (seed,seed))
    file.close()
