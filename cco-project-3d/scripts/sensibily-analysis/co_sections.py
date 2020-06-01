Q_PERF = 8.33e-06
P_PERF = 1.33e+04
P_TERM = 9.58e+03
V_PERF = 1.0e-04
N_TERM = 130
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

def write_co_main_section (file,seed,rand_offset,start_radius):
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
    file.write("use_only_murray = true\n")
    file.write("start_radius = %g\n" % start_radius)
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

def write_co_custom_cost_function_section (file):
    file.write("[cost_function]\n")
    file.write("library_name = shared-libs/libminimize_custom_function.so\n")
    file.write("function_name = minimize_custom_function\n")
    file.write("beta = 1.0\n")
    file.write("alpha = 0.0\n")
    file.write("\n")

def write_co_activation_time_cost_function_section (file):
    file.write("[cost_function]\n")
    file.write("library_name = shared-libs/libminimize_activation_time.so\n")
    file.write("function_name = minimize_tree_activation_time\n")
    file.write("G = 7.9\n")
    file.write("Cf = 3.4\n")
    file.write("tauf = 0.1\n")
    file.write("\n")

def write_cco_pruning_section (file):
    file.write("[pruning]\n")
    file.write("use_pruning = true\n")
    file.write("A = %g\n" % A)
    file.write("B = %g\n" % B)
    file.write("C = %g\n" % C)
    file.write("\n")

def write_co_config_file (seed,rand_offset):
    filename = "cco_config/elizabeth_purkinje_co_seed-%u_offset-%u_nterm-130.ini" % (seed,rand_offset)
    file = open(filename,"w")

    write_cco_main_section(file,seed,rand_offset)    
    write_cco_cloud_points_section(file)
    write_cco_local_optimization_section(file)
    write_co_custom_cost_function_section(file)
    write_cco_pruning_section(file)

    file.close()

def write_co_activation_time_config_file (seed,rand_offset,start_radius):
    filename = "cco_config/elizabeth_purkinje_co_activation_time_seed-%u_offset-%u_nterm-130.ini" % (seed,rand_offset)
    file = open(filename,"w")

    write_co_main_section(file,seed,rand_offset,start_radius)    
    write_cco_cloud_points_section(file)
    write_cco_local_optimization_section(file)
    write_co_activation_time_cost_function_section(file)
    write_cco_pruning_section(file)

    file.close()
