# =====================================================================================================
# This library configures the solver to generate Purkinje network that follows the cost function:
#   MINIMIZE WEIGTHED ACTIVATION TIME WITH ANGLE RESTRICTION 
# =====================================================================================================

Q_PERF = 8.33e-06
P_PERF = 1.33e+04
P_TERM = 9.58e+03
V_PERF = 1.0e-04
N_TERM = 80
ROOT_X = 0.0010
ROOT_Y = 0.0985
ROOT_Z = 0.0005
START_RADIUS = 0.1
GAMMA = 19.0

SAVE_NETWORK_PATH = "outputs/02_Oxford/04_Minimize_Total_Linked"

USE_CLOUD_POINTS = True
CLOUD_POINTS_FILENAME = "clouds/generated-clouds/slab_coordinates_guided.pts"
USE_OBSTACLE = False
USE_PMJ_LOCATION = True
PMJ_LOCATION_FILENAME = "clouds/generated-clouds/slab_pmj_locations.pts"

USE_LOCAL_OPTIMIZATION = True
LOCAL_OPTIMIZATION_FUNCTION = "rafael_local_optimization"

COST_FUNCTION_LIBRARY = "shared-libs/libminimize_activation_time.so"
COST_FUNCTION_NAME = "minimize_weighted_activation_time_with_angle_restriction"
W = 0.0
MIN_DEGREES_LIMIT = 1.0
MAX_DEGREES_LIMIT = 63.0

USE_PRUNING = False
PRUNING_FUNCTION_NAME = "hyperbolic_tangent"
PRUNING_PARAM_A = 50.0
PRUNING_PARAM_B = -0.25
PRUNING_PARAM_C = 3.0
PRUNING_PARAM_D = 50.0
#PRUNING_FUNCTION_NAME = "length"
#PRUNING_PARAM_LENGTH_LIMIT = 0.01

def write_co_main_section (file,seed,rand_offset):
    file.write("[main]\n")
    file.write("Q_perf = %g\n" % Q_PERF)
    file.write("p_perf = %g\n" % P_PERF)
    file.write("p_term = %g\n" % P_TERM)
    file.write("V_perf = %g\n" % V_PERF)
    file.write("N_term = %u\n" % N_TERM)
    file.write("root_x = %g\n" % ROOT_X)
    file.write("root_y = %g\n" % ROOT_Y)
    file.write("root_z = %g\n" % ROOT_Z)
    file.write("max_rand_offset = %u\n" % rand_offset)
    file.write("seed = %u\n" % seed)
    file.write("use_only_murray = true\n")
    file.write("start_radius = %g\n" % START_RADIUS)
    file.write("gamma = %g\n" % GAMMA)
    file.write("\n")

def write_co_save_network (file,seed):
    file.write("[save_network]\n")
    file.write("output_dir = %s/seed:%d\n" % (SAVE_NETWORK_PATH,seed))
    file.write("\n")

def write_co_cloud_points_section (file):
    file.write("[cloud_points]\n")
    file.write("use_cloud_points = true\n")
    file.write("cloud_points_filename = %s\n" % CLOUD_POINTS_FILENAME)
    file.write("use_obstacle = false\n")
    if (USE_PMJ_LOCATION):
        file.write("use_pmj_location = true\n")
        file.write("pmj_location_filename = %s\n" % (PMJ_LOCATION_FILENAME))
    else:
        file.write("use_pmj_location = false\n")
    file.write("\n")

def write_co_local_optimization_section (file):
    file.write("[local_optimization]\n")
    file.write("use_local_optimization = true\n")
    file.write("local_optimization_function = %s\n" % LOCAL_OPTIMIZATION_FUNCTION)
    file.write("\n")

def write_co_custom_cost_function_section (file):
    file.write("[cost_function]\n")
    file.write("library_name = %s\n" % (COST_FUNCTION_LIBRARY)) 
    file.write("function_name = %s\n" % (COST_FUNCTION_NAME))
    file.write("w = %g\n" % (W))
    file.write("G = 7.9\n")
    file.write("Cf = 3.4\n")
    file.write("tauf = 0.1\n")
    file.write("min_degrees_limit = %g\n" % (MIN_DEGREES_LIMIT))
    file.write("max_degrees_limit = %g\n" % (MAX_DEGREES_LIMIT))
    file.write("\n")

def write_co_pruning_section (file):
    file.write("[pruning]\n")
    file.write("use_pruning = true\n")
    file.write("pruning_function = %s\n" % (PRUNING_FUNCTION_NAME))
    #file.write("length_limit = %g\n" % (PRUNING_PARAM_LENGTH_LIMIT))
    file.write("A = %s\n" % (PRUNING_PARAM_A))
    file.write("B = %s\n" % (PRUNING_PARAM_B))
    file.write("C = %s\n" % (PRUNING_PARAM_C))
    file.write("D = %s\n" % (PRUNING_PARAM_D))
    file.write("\n")

def write_co_activation_time_config_file (seed,rand_offset):

    filename = "network/min:total-linked_seed-%u.ini" % (seed)
    file = open(filename,"w")

    write_co_main_section(file,seed,rand_offset)    
    write_co_save_network(file,seed)
    write_co_cloud_points_section(file)
    write_co_local_optimization_section(file)
    write_co_custom_cost_function_section(file)
    if (USE_PRUNING):
        write_co_pruning_section(file)

    file.close()
