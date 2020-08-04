# =====================================================================================================
# This library configures the solver to generate Purkinje network that follows the cost function:
#   MINIMIZE TOTAL LENGTH WITH ANGLE RESTRICTION (the same one used on the IEEE paper from 2018)
# =====================================================================================================

Q_PERF = 8.33e-06
P_PERF = 1.33e+04
P_TERM = 9.58e+03
V_PERF = 1.0e-04
N_TERM = 650
ROOT_X = -0.004
ROOT_Y = 0.02
ROOT_Z = 0.0265
START_RADIUS = 0.00102269

SAVE_NETWORK_PATH = "outputs/01_SRN_Purkinje/04_CO_Length_With_Pruning"
# SAVE_NETWORK_PATH = "outputs/01_SRN_Purkinje/01_CO_Length" 

USE_CLOUD_POINTS = True
CLOUD_POINTS_FILENAME = "clouds/private/elizabeth_remapped_guided_3.pts"
USE_OBSTACLE = False
USE_PMJ_LOCATION = False

USE_LOCAL_OPTIMIZATION = True
LOCAL_OPTIMIZATION_FUNCTION = "rafael_local_optimization"

COST_FUNCTION_LIBRARY = "shared-libs/libminimize_custom_function.so"
COST_FUNCTION_NAME = "minimize_custom_function_with_angle_restriction"
MIN_DEGREES_LIMIT = 1.0
MAX_DEGREES_LIMIT = 63.0

USE_PRUNING = True
PRUNING_FUNCTION_NAME = "hyperbolic_tangent"
PRUNING_PARAM_A = 50.0
PRUNING_PARAM_B = -0.25
PRUNING_PARAM_C = 3.0
PRUNING_PARAM_D = 50.0

def write_co_main_section (file,seed,rand_offset):
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

def write_co_save_network (file,seed):
    file.write("[save_network]\n")
    file.write("output_dir = %s/seed:%d\n" % (SAVE_NETWORK_PATH,seed))
    file.write("\n")

def write_co_cloud_points_section (file):
    file.write("[cloud_points]\n")
    file.write("use_cloud_points = true\n")
    file.write("cloud_points_filename = %s\n" % CLOUD_POINTS_FILENAME)
    file.write("use_obstacle = false\n")
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
    file.write("beta = 1.0\n")
    file.write("alpha = 0.0\n")
    file.write("min_degrees_limit = %g\n" % (MIN_DEGREES_LIMIT))
    file.write("max_degrees_limit = %g\n" % (MAX_DEGREES_LIMIT))
    file.write("\n")

def write_co_pruning_section (file):
    file.write("[pruning]\n")
    file.write("use_pruning = true\n")
    file.write("pruning_function = %s\n" % (PRUNING_FUNCTION_NAME))
    file.write("A = %s\n" % (PRUNING_PARAM_A))
    file.write("B = %s\n" % (PRUNING_PARAM_B))
    file.write("C = %s\n" % (PRUNING_PARAM_C))
    file.write("D = %s\n" % (PRUNING_PARAM_D))
    file.write("\n")

def write_co_config_file (seed,rand_offset):
    filename = "network/co_min:length_with_pruning_seed-%u_nterm:650.ini" % (seed)
    file = open(filename,"w")

    write_co_main_section(file,seed,rand_offset)    
    write_co_save_network(file,seed)
    write_co_cloud_points_section(file)
    write_co_local_optimization_section(file)
    write_co_custom_cost_function_section(file)
    if (USE_PRUNING):
        write_co_pruning_section(file)

    file.close()
