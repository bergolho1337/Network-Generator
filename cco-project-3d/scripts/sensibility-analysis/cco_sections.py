Q_PERF = 8.33e-06
P_PERF = 1.33e+04
P_TERM = 9.58e+03
V_PERF = 1.0e-04
N_TERM = 200
ROOT_X = -0.004
ROOT_Y = 0.02
ROOT_Z = 0.0265
START_RADIUS = 0.00102269
A = 100.0
B = -0.25
C = 0.0

USE_CLOUD_POINTS = True
CLOUD_POINTS_FILENAME = "clouds/private/elizabeth_exterior_LV_remapped.pts"
USE_OBSTACLE = False
OBSTACLE_FILENAME = "clouds/private/elizabeth_interior_LV.stl"
PMJ_LOCATION_FILENAME = "clouds/private/elizabeth-pmj-location_LV.pts"

USE_LOCAL_OPTIMIZATION = True
LOCAL_OPTIMIZATION_FUNCTION = "rafael_local_optimization"

COST_FUNCTION_LIBRARY_NAME = "shared-libs/libminimize_volume.so"
COST_FUNCTION_NAME = "minimize_tree_volume_default"

def write_cco_main_section (file,seed,rand_offset):
    file.write("[main]\n")
    file.write("Q_perf = %g\n" % Q_PERF)
    file.write("p_perf = %g\n" % P_PERF)
    file.write("p_term = %g\n" % P_TERM)
    file.write("p_perf = %g\n" % P_PERF)
    file.write("V_perf = %g\n" % V_PERF)
    file.write("N_term = %u\n" % N_TERM)
    file.write("root_x = %g\n" % ROOT_X)
    file.write("root_y = %g\n" % ROOT_Y)
    file.write("root_z = %g\n" % ROOT_Z)
    file.write("max_rand_offset = %u\n" % rand_offset)
    file.write("seed = %u\n" % seed)
    file.write("use_only_murray = false\n")
    file.write("\n")

def write_cco_save_network_section(file,seed,rand_offset):
    file.write("[save_network]\n")
    file.write("output_dir = outputs/elizabeth_minimize_volume_seed:%u_offset:%u\n" % (seed,rand_offset))
    file.write("\n")

def write_cco_cloud_points_section (file):
    file.write("[cloud_points]\n")
    file.write("use_cloud_points = true\n")
    file.write("cloud_points_filename = %s\n" % CLOUD_POINTS_FILENAME)
    file.write("use_obstacle = false\n")
    file.write("obstacle_filename = %s\n" % OBSTACLE_FILENAME)
    file.write("use_pmj_location = true\n")
    file.write("pmj_location_filename = %s\n" % PMJ_LOCATION_FILENAME)
    file.write("\n")

def write_cco_local_optimization_section (file):
    file.write("[local_optimization]\n")
    file.write("use_local_optimization = true\n")
    file.write("local_optimization_function = %s\n" % LOCAL_OPTIMIZATION_FUNCTION)
    file.write("\n")

def write_cco_cost_function_section (file):
    file.write("[cost_function]\n")
    file.write("library_name = %s\n" % COST_FUNCTION_LIBRARY_NAME)
    file.write("function_name = %s\n" % COST_FUNCTION_NAME)
    file.write("\n")

def write_cco_pruning_section (file):
    file.write("[pruning]\n")
    file.write("use_pruning = true\n")
    file.write("A = %g\n" % A)
    file.write("B = %g\n" % B)
    file.write("C = %g\n" % C)
    file.write("\n")

def write_cco_config_file (seed,rand_offset):
    filename = "minimize_volume/elizabeth_seed-%u_offset-%u.ini" % (seed,rand_offset)
    file = open(filename,"w")

    write_cco_main_section(file,seed,rand_offset)
    write_cco_save_network_section(file,seed,rand_offset)    
    write_cco_cloud_points_section(file)
    write_cco_local_optimization_section(file)
    write_cco_cost_function_section(file)
    write_cco_pruning_section(file)

    file.close()