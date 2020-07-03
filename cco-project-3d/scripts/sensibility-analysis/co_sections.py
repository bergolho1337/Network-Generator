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
A = 100.0
B = -0.25
C = 0.0

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
    file.write("output_dir = outputs/01_SRN/01_CO_Length/seed:%d\n" % (seed))
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

'''
def write_co_activation_time_cost_function_section (file):
    file.write("[cost_function]\n")
    file.write("library_name = shared-libs/libminimize_activation_time.so\n")
    file.write("function_name = minimize_tree_activation_time\n")
    file.write("G = 7.9\n")
    file.write("Cf = 3.4\n")
    file.write("tauf = 0.1\n")
    file.write("\n")
'''

def write_co_pruning_section (file):
    file.write("[pruning]\n")
    file.write("use_pruning = false\n")
    file.write("\n")

def write_co_config_file (seed,rand_offset):
    filename = "ieee/co_min:length_seed-%u_nterm:650.ini" % (seed)
    file = open(filename,"w")

    write_co_main_section(file,seed,rand_offset)    
    write_co_save_network(file,seed)
    write_co_cloud_points_section(file)
    write_co_local_optimization_section(file)
    write_co_custom_cost_function_section(file)
    write_co_pruning_section(file)

    file.close()

'''
def write_co_activation_time_config_file (seed,rand_offset,start_radius):
    filename = "cco_config/elizabeth_purkinje_co_activation_time_seed-%u_offset-%u_nterm-130.ini" % (seed,rand_offset)
    file = open(filename,"w")

    write_co_main_section(file,seed,rand_offset,start_radius)    
    write_cco_cloud_points_section(file)
    write_cco_local_optimization_section(file)
    write_co_activation_time_cost_function_section(file)
    write_cco_pruning_section(file)

    file.close()
'''