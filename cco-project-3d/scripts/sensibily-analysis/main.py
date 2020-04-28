import sys
import chaospy as cp

Q_PERF = 8.33e-06
P_PERF = 1.33e+04
P_TERM = 9.58e+03
V_PERF = 1.0e-04
N_TERM = 130
ROOT_X = -0.004
ROOT_Y = 0.02
ROOT_Z = 0.0265
USE_ONLY_MURRAY = False

USE_CLOUD_POINTS = True
CLOUD_POINTS_FILENAME = "clouds/private/elizabeth_exterior_LV_remapped.pts"
USE_OBSTACLE = False
OBSTACLE_FILENAME = "clouds/private/elizabeth_interior_LV.stl"

USE_LOCAL_OPTIMIZATION = True
LOCAL_OPTIMIZATION_FUNCTION = "rafael_local_optimization"

COST_FUNCTION_LIBRARY_NAME = "shared-libs/libminimize_volume.so"
COST_FUNCTION_NAME = "minimize_tree_volume_default"

def write_cco_config_file (seed,rand_offset):
    filename = "cco_config/elizabeth_purkinje_cco_seed-%u_offset-%u_nterm-130.ini" % (seed,rand_offset)
    file = open(filename,"w")

    # Write [main] section
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

    # Write [cloud_points] section
    file.write("[cloud_points]\n")
    file.write("use_cloud_points = true\n")
    file.write("cloud_points_filename = %s\n" % CLOUD_POINTS_FILENAME)
    file.write("use_obstacle = false\n")
    file.write("obstacle_filename = %s\n" % OBSTACLE_FILENAME)
    file.write("\n")

    # Write [local_optimization] section
    file.write("[local_optimization]\n")
    file.write("use_local_optimization = true\n")
    file.write("local_optimization_function = %s\n" % LOCAL_OPTIMIZATION_FUNCTION)
    file.write("\n")

    # Write [cost_function] section
    file.write("[cost_function]\n")
    file.write("library_name = %s\n" % COST_FUNCTION_LIBRARY_NAME)
    file.write("function_name = %s\n" % COST_FUNCTION_NAME)
    file.write("\n")

    file.close()

def write_monoalg_config_file (seed,rand_offset):
    filename = "monoalg3d_config/elizabeth_purkinje_cco_seed-%u_offset-%u_nterm-130.ini" % (seed,rand_offset)
    file = open(filename,"w")

    # Write [main] section
    file.write("[main]\n")
    file.write("num_threads = 4\n")
    file.write("dt_pde = 0.02\n")
    file.write("simulation_time = 800.0\n")
    file.write("abort_on_no_activity = false\n")
    file.write("use_adaptivity = false\n")
    file.write("calc_activation_time = true\n")
    file.write("print_conductivity_map = false\n")
    file.write("\n")

    # Write [update_monodomain] section
    file.write("[update_monodomain]\n")
    file.write("main_function = update_monodomain_default\n")
    #file.write("library_file=/home/berg/MonoAlg3D_C/shared_libs/libdefault_update_monodomain.so\n")
    file.write("\n")

    # Write [save_result] section
    file.write("[save_result]\n")
    file.write("print_rate = 100\n")
    file.write("output_dir = ./outputs/elizabeth_coupled_arpf_cco_seed-%u_offset-%u_LV\n" % (seed,rand_offset))
    file.write("main_function = save_as_vtu_tissue_coupled_vtp_purkinje\n")
    file.write("save_pvd = true\n")
    file.write("file_prefix = V_Tissue\n")
    file.write("file_prefix_purkinje = V_Purkinje\n")
    file.write("binary = false\n")
    file.write("compress = false\n")
    #file.write("library_file=/home/berg/MonoAlg3D_C/shared_libs/libdefault_save_mesh.so\n")
    file.write("\n")

    # Write [assembly_matrix] section
    file.write("[assembly_matrix]\n")
    file.write("init_function = set_initial_conditions_coupled_fvm\n")
    file.write("sigma_x = 0.00005336\n")
    file.write("sigma_y = 0.00005336\n")
    file.write("sigma_z = 0.00005336\n")
    file.write("sigma_purkinje = 0.004\n")
    file.write("main_function = purkinje_coupled_endocardium_assembly_matrix\n")
    file.write("library_file = shared_libs/libpurkinje_coupled_matrix_assembly.so\n")
    #file.write("library_file=/home/berg/MonoAlg3D_C/shared_libs/libpurkinje_coupled_matrix_assembly.so\n")
    file.write("\n")

    # Write [linear_system_solver] section
    file.write("[linear_system_solver]\n")
    file.write("tolerance = 1e-16\n")
    file.write("use_preconditioner = yes\n")
    file.write("max_iterations = 200\n")
    file.write("main_function = conjugate_gradient\n")
    file.write("library_file = shared_libs/libdefault_linear_system_solver.so\n")
    #file.write("library_file = /home/berg/MonoAlg3D_C/shared_libs/libdefault_linear_system_solver.so\n")
    file.write("\n")

    # Write [alg] section
    file.write("[alg]\n")
    file.write("refinement_bound = 0.11\n")
    file.write("derefinement_bound = 0.10\n")
    file.write("refine_each = 1\n")
    file.write("derefine_each = 1\n")
    file.write("\n")

    # Write [domain] section
    file.write("[domain]\n")
    file.write("name = Elizabeth LV Canine Endocardium Mesh\n")
    file.write("maximum_discretization = 500.0\n")
    file.write("start_dx = 500.0\n")
    file.write("start_dy = 500.0\n")
    file.write("start_dz = 500.0\n")
    file.write("x_domain_limit = 128000.0\n")
    file.write("y_domain_limit = 128000.0\n")
    file.write("z_domain_limit = 128000.0\n")
    file.write("refinement_steps = 7\n")
    file.write("total_number_mesh_points = 1061776\n")
    file.write("mesh_file = meshes/elizabeth-canine-lv-endocardium.alg\n")
    file.write("main_function = initialize_grid_with_custom_mesh\n")
    #file.write("library_file = /home/berg/MonoAlg3D_C/shared_libs/libdefault_domains.so.so\n")    
    file.write("\n")

    # Write [purkinje] section
    file.write("[purkinje]\n")
    file.write("name = Simple Purkinje\n")
    file.write("rpmj = 1.0e+02\n")
    file.write("pmj_scale = 0.01\n")
    file.write("start_discretization = 100.0\n")
    file.write("start_dx = 100.0\n")
    file.write("retro_propagation = false\n")
    file.write("network_file = networks/elizabeth-meshes/cco_classic/elizabeth_purkinje_cco_seed-%u_offset-%u_nterm-130.vtk\n" % (seed,rand_offset))
    file.write("main_function = initialize_purkinje_with_custom_mesh\n")
    file.write("library_file = shared_libs/libdefault_purkinje.so\n")
    #file.write("library_file = /home/berg/MonoAlg3D_C/shared_libs/libdefault_purkinje.so\n")
    file.write("\n")

    # Write [purkinje_ode_solver] section
    file.write("[purkinje_ode_solver]\n")
    file.write("dt_ode = 0.01\n")
    file.write("use_gpu = no\n")
    file.write("gpu_id = 0\n")
    file.write("library_file = shared_libs/libstewart_aslanidi_noble_2009.so\n")
    #file.write("library_file = /home/berg/MonoAlg3D_C/shared_libs/libstewart_aslanidi_noble_2009.so\n")
    file.write("\n")

    # Write [ode_solver] section
    file.write("[ode_solver]\n")
    file.write("dt_ode = 0.02\n")
    file.write("use_gpu = yes\n")
    file.write("gpu_id = 0\n")
    file.write("library_file = shared_libs/libten_tusscher_2006.so\n")
    #file.write("library_file = /home/berg/MonoAlg3D_C/shared_libs/libten_tusscher_2006.so\n")
    file.write("\n")

    # Write [stimulus] section
    file.write("[stim_purkinje_his]\n")
    file.write("start = 0.0\n")
    file.write("duration = 4.0\n")
    file.write("current = -50.0\n")
    file.write("id_limit = 20\n")
    file.write("main_function = stim_if_id_less_than\n")
    #file.write("library_file = /home/berg/MonoAlg3D_C/shared_libs/libdefault_stimuli.so")
    file.write("\n")

    file.close()

def main():
	
    seeds = [1562046115,1562013988,1562042299,1562005512,1562009134,1562009768,1562044566,1562008423,1562036996,1562020974]
    rand_offsets = [2,3,4,5,6,7,8,9,10]

    for seed in seeds:
        for rand_offset in rand_offsets:
            write_cco_config_file(seed,rand_offset)
            write_monoalg_config_file(seed,rand_offset)

if __name__ == "__main__":
	main()


'''
# ChaosPy code

min_rand_offset = 1
max_rand_offset = 10

rand_offset = cp.Uniform(min_rand_offset,max_rand_offset)
#b = cp.Uniform(0,5)

#distribution = cp.J(a,b)
distribution = cp.J(rand_offset)

samples = distribution.sample(10,"L")

print samples.T.astype(int)
'''
