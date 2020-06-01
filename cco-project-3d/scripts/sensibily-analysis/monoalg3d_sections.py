def write_monoalg_main_section (file,seed,rand_offset):
    file.write("[main]\n")
    file.write("num_threads = 6\n")
    file.write("dt_pde = 0.02\n")
    file.write("simulation_time = 500.0\n")
    file.write("abort_on_no_activity = false\n")
    file.write("use_adaptivity = false\n")
    file.write("calc_activation_time = true\n")
    file.write("print_conductivity_map = false\n")
    file.write("\n")

def write_monoalg_update_monodomain (file,use_cluster=False):
    file.write("[update_monodomain]\n")
    file.write("main_function = update_monodomain_default\n")
    if (use_cluster):
        file.write("library_file=/home/berg/MonoAlg3D_C/shared_libs/libdefault_update_monodomain.so\n")
    file.write("\n")

def write_monoalg_save_result_section (file,seed,rand_offset,algorithm_number,use_cluster=False):
    file.write("[save_result]\n")
    file.write("print_rate = 1000\n")
    if (algorithm_number == 1):
        file.write("output_dir = ./outputs/elizabeth_coupled_arpf_cco_seed-%u_offset-%u_LV\n" % (seed,rand_offset))
    elif (algorithm_number == 2):
        file.write("output_dir = ./outputs/elizabeth_coupled_arpf_co_seed-%u_offset-%u_LV\n" % (seed,rand_offset))
    elif (algorithm_number == 3):
        file.write("output_dir = ./outputs/elizabeth_coupled_arpf_co_activation_time_seed-%u_offset-%u_LV\n" % (seed,rand_offset))
    file.write("main_function = save_as_vtu_tissue_coupled_vtp_purkinje\n")
    file.write("save_pvd = true\n")
    file.write("file_prefix = V_Tissue\n")
    file.write("file_prefix_purkinje = V_Purkinje\n")
    file.write("binary = false\n")
    file.write("compress = false\n")
    if (use_cluster):
        file.write("library_file=/home/berg/MonoAlg3D_C/shared_libs/libdefault_save_mesh.so\n")
    file.write("\n")

def write_monoalg_assembly_matrix_section (file,use_cluster=False):
    file.write("[assembly_matrix]\n")
    file.write("init_function = set_initial_conditions_coupled_fvm\n")
    file.write("sigma_x = 0.00005336\n")
    file.write("sigma_y = 0.00005336\n")
    file.write("sigma_z = 0.00005336\n")
    file.write("sigma_purkinje = 0.004\n")
    file.write("main_function = purkinje_coupled_endocardium_assembly_matrix\n")
    file.write("library_file = shared_libs/libpurkinje_coupled_matrix_assembly.so\n")
    if (use_cluster):
        file.write("library_file=/home/berg/MonoAlg3D_C/shared_libs/libpurkinje_coupled_matrix_assembly.so\n")
    file.write("\n")

def write_monoalg_linear_system_solver_section (file,use_cluster=False):
    file.write("[linear_system_solver]\n")
    file.write("tolerance = 1e-16\n")
    file.write("use_preconditioner = no\n")
    file.write("max_iterations = 200\n")
    file.write("main_function = conjugate_gradient\n")
    file.write("library_file = shared_libs/libdefault_linear_system_solver.so\n")
    if (use_cluster):
        file.write("library_file = /home/berg/MonoAlg3D_C/shared_libs/libdefault_linear_system_solver.so\n")
    file.write("\n")

def write_monoalg_alg_section (file,use_cluster=False):
    file.write("[alg]\n")
    file.write("refinement_bound = 0.11\n")
    file.write("derefinement_bound = 0.10\n")
    file.write("refine_each = 1\n")
    file.write("derefine_each = 1\n")
    file.write("\n")

def write_monoalg_domain_section (file,use_cluster=False):
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
    if (use_cluster):
        file.write("library_file = /home/berg/MonoAlg3D_C/shared_libs/libdefault_domains.so.so\n")    
    file.write("\n")

def write_monoalg_purkinje_section (file,seed,rand_offset,algorithm_number,use_cluster=False):
    file.write("[purkinje]\n")
    file.write("name = Simple Purkinje\n")
    file.write("rpmj = 1.0e+02\n")
    file.write("pmj_scale = 0.01\n")
    file.write("start_discretization = 200.0\n")
    file.write("start_dx = 200.0\n")
    file.write("retro_propagation = true\n")
    if (algorithm_number == 1):
        file.write("network_file = networks/elizabeth-meshes/cco_classic/elizabeth_purkinje_cco_seed-%u_offset-%u_nterm-130.vtk\n" % (seed,rand_offset))
    elif (algorithm_number == 2):
        file.write("network_file = networks/elizabeth-meshes/co_length/elizabeth_purkinje_co_seed-%u_offset-%u_nterm-130.vtk\n" % (seed,rand_offset))
    elif (algorithm_number == 3):
        file.write("network_file = networks/elizabeth-meshes/co_activation_time/elizabeth_purkinje_co_activation_time_seed-%u_offset-%u_nterm-130.vtk\n" % (seed,rand_offset))
    file.write("main_function = initialize_purkinje_with_custom_mesh\n")
    file.write("library_file = shared_libs/libdefault_purkinje.so\n")
    if (use_cluster):
        file.write("library_file = /home/berg/MonoAlg3D_C/shared_libs/libdefault_purkinje.so\n")
    file.write("\n")

def write_monoalg_purkinje_ode_solver_section (file,use_cluster=False):
    file.write("[purkinje_ode_solver]\n")
    file.write("dt_ode = 0.02\n")
    file.write("use_gpu = no\n")
    file.write("gpu_id = 0\n")
    file.write("library_file = shared_libs/libstewart_aslanidi_noble_2009.so\n")
    if (use_cluster):
        file.write("library_file = /home/berg/MonoAlg3D_C/shared_libs/libstewart_aslanidi_noble_2009.so\n")
    file.write("\n")

def write_monoalg_ode_solver_section (file,use_cluster=False):
    file.write("[ode_solver]\n")
    file.write("dt_ode = 0.02\n")
    file.write("use_gpu = yes\n")
    file.write("gpu_id = 0\n")
    file.write("library_file = shared_libs/libten_tusscher_2006.so\n")
    if (use_cluster):
        file.write("library_file = /home/berg/MonoAlg3D_C/shared_libs/libten_tusscher_2006.so\n")
    file.write("\n")

def write_monoalg_stimulus_section (file,use_cluster=False):
    file.write("[stim_purkinje_his]\n")
    file.write("start = 0.0\n")
    file.write("duration = 2.0\n")
    file.write("current = -90.0\n")
    file.write("id_limit = 5\n")
    file.write("main_function = stim_if_id_less_than\n")
    if (use_cluster):
        file.write("library_file = /home/berg/MonoAlg3D_C/shared_libs/libdefault_stimuli.so")
    file.write("\n")

'''
Network tester
'''
'''
def write_monoalg_main_section (file,seed,rand_offset):
    file.write("[main]\n")
    file.write("num_threads = 6\n")
    file.write("dt_pde = 0.02\n")
    file.write("simulation_time = 500.0\n")
    file.write("abort_on_no_activity = false\n")
    file.write("use_adaptivity = false\n")
    file.write("calc_activation_time = true\n")
    file.write("print_conductivity_map = false\n")
    file.write("\n")

def write_monoalg_update_monodomain (file,use_cluster=False):
    file.write("[update_monodomain]\n")
    file.write("main_function = update_monodomain_default\n")
    if (use_cluster):
        file.write("library_file=/home/berg/MonoAlg3D_C/shared_libs/libdefault_update_monodomain.so\n")
    file.write("\n")

def write_monoalg_save_result_section (file,seed,rand_offset,use_cluster=False):
    file.write("[save_result]\n")
    file.write("print_rate = 1000\n")
    file.write("output_dir = ./outputs/elizabeth_coupled_arpf_cco_seed-%u_offset-%u_LV\n" % (seed,rand_offset))
    file.write("main_function = save_as_vtp_purkinje\n")
    file.write("save_pvd = true\n")
    file.write("file_prefix = V\n")
    file.write("binary = false\n")
    file.write("compress = false\n")
    if (use_cluster):
        file.write("library_file=/home/berg/MonoAlg3D_C/shared_libs/libdefault_save_mesh.so\n")
    file.write("\n")

def write_monoalg_assembly_matrix_section (file,use_cluster=False):
    file.write("[assembly_matrix]\n")
    file.write("init_function = set_initial_conditions_fvm\n")
    file.write("sigma_x = 0.004\n")
    file.write("sigma_y = 0.004\n")
    file.write("sigma_z = 0.004\n")
    file.write("sigma_purkinje = 0.004\n")
    file.write("main_function = purkinje_fibers_assembly_matrix\n")
    file.write("library_file = shared_libs/libpurkinje_matrix_assembly.so\n")
    if (use_cluster):
        file.write("library_file=/home/berg/MonoAlg3D_C/shared_libs/libpurkinje_coupled_matrix_assembly.so\n")
    file.write("\n")

def write_monoalg_linear_system_solver_section (file,use_cluster=False):
    file.write("[linear_system_solver]\n")
    file.write("tolerance = 1e-16\n")
    file.write("use_preconditioner = no\n")
    file.write("max_iterations = 200\n")
    file.write("main_function = conjugate_gradient\n")
    file.write("library_file = shared_libs/libdefault_linear_system_solver.so\n")
    if (use_cluster):
        file.write("library_file = /home/berg/MonoAlg3D_C/shared_libs/libdefault_linear_system_solver.so\n")
    file.write("\n")

def write_monoalg_purkinje_section (file,seed,rand_offset,use_cluster=False):
    file.write("[purkinje]\n")
    file.write("name = Simple Purkinje\n")
    file.write("start_discretization = 200.0\n")
    file.write("start_dx = 200.0\n")
    file.write("network_file = networks/elizabeth-meshes/cco_classic/elizabeth_purkinje_cco_seed-%u_offset-%u_nterm-130.vtk\n" % (seed,rand_offset))
    file.write("main_function = initialize_purkinje_with_custom_mesh\n")
    file.write("library_file = shared_libs/libdefault_purkinje.so\n")
    if (use_cluster):
        file.write("library_file = /home/berg/MonoAlg3D_C/shared_libs/libdefault_purkinje.so\n")
    file.write("\n")

def write_monoalg_ode_solver_section (file,use_cluster=False):
    file.write("[ode_solver]\n")
    file.write("dt_ode = 0.02\n")
    file.write("use_gpu = no\n")
    file.write("gpu_id = 0\n")
    file.write("library_file = shared_libs/libstewart_aslanidi_noble_2009.so\n")
    if (use_cluster):
        file.write("library_file = /home/berg/MonoAlg3D_C/shared_libs/libten_tusscher_2006.so\n")
    file.write("\n")

def write_monoalg_stimulus_section (file,use_cluster=False):
    file.write("[stim_purkinje_his]\n")
    file.write("start = 0.0\n")
    file.write("duration = 2.0\n")
    file.write("current = -90.0\n")
    file.write("id_limit = 5\n")
    file.write("main_function = stim_if_id_less_than\n")
    if (use_cluster):
        file.write("library_file = /home/berg/MonoAlg3D_C/shared_libs/libdefault_stimuli.so")
    file.write("\n")
'''

def write_monoalg_config_file (seed,rand_offset,algorithm_number):
    filename = "monoalg3d_config/elizabeth_purkinje_cco_seed-%u_offset-%u_nterm-130.ini" % (seed,rand_offset)
    file = open(filename,"w")

    write_monoalg_main_section(file,seed,rand_offset)
    write_monoalg_update_monodomain(file)
    write_monoalg_save_result_section(file,seed,rand_offset,algorithm_number)
    write_monoalg_assembly_matrix_section(file)
    write_monoalg_linear_system_solver_section(file)
    write_monoalg_alg_section(file)
    write_monoalg_domain_section(file)
    write_monoalg_purkinje_section(file,seed,rand_offset,algorithm_number)
    write_monoalg_purkinje_ode_solver_section(file)
    write_monoalg_ode_solver_section(file)
    write_monoalg_stimulus_section(file)

    file.close()
