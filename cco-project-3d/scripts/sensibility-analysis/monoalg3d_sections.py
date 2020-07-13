def write_monoalg_main_section (file,seed,rand_offset):
    file.write("[main]\n")
    file.write("num_threads = 10\n")
    file.write("dt_pde = 0.02\n")
    file.write("simulation_time = 500.0\n")
    file.write("abort_on_no_activity = false\n")
    file.write("use_adaptivity = false\n")
    file.write("\n")

def write_monoalg_update_monodomain (file,use_cluster=False):
    file.write("[update_monodomain]\n")
    file.write("main_function = update_monodomain_default\n")
    if (use_cluster):
        file.write("library_file=/home/berg/MonoAlg3D_C/shared_libs/libdefault_update_monodomain.so\n")
    file.write("\n")

def write_monoalg_save_result_section (file,seed,rand_offset,algorithm_number,use_cluster=False):
    file.write("[save_result]\n")
    file.write("print_rate = 25\n")
    #file.write("output_dir = ./outputs/elizabeth_co:min_length_seed:%u\n" % (seed))
    file.write("output_dir = ./outputs/elizabeth_co:min_at_seed:%u\n" % (seed))
    file.write("init_function=init_save_purkinje_coupling_with_activation_times\n")
    file.write("end_function=end_save_purkinje_coupling_with_activation_times\n")
    file.write("main_function=save_purkinje_coupling_with_activation_times\n")
    file.write("apd_threshold_tissue=-70.0\n")
    file.write("apd_threshold_purkinje=-70.0\n")
    file.write("save_activation_time=true\n")
    file.write("save_apd=true\n")
    file.write("save_pvd = true\n")
    file.write("file_prefix = V_Tissue\n")
    file.write("file_prefix_purkinje = V_Purkinje\n")
    file.write("remove_older_simulation=true\n")
    if (use_cluster):
        file.write("library_file=/home/berg/MonoAlg3D_C/shared_libs/libdefault_save_mesh.so\n")
    file.write("\n")

def write_monoalg_assembly_matrix_section (file,use_cluster=False):
    file.write("[assembly_matrix]\n")
    file.write("init_function = set_initial_conditions_coupling_fvm\n")
    file.write("sigma_x = 0.00005336\n")
    file.write("sigma_y = 0.00005336\n")
    file.write("sigma_z = 0.00005336\n")
    file.write("sigma_purkinje = 0.004\n")
    file.write("main_function = purkinje_coupling_assembly_matrix\n")
    file.write("library_file = shared_libs/libpurkinje_coupling_matrix_assembly.so\n")
    if (use_cluster):
        file.write("library_file=/home/berg/MonoAlg3D_C/shared_libs/libpurkinje_coupling_matrix_assembly.so\n")
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
    file.write("mesh_file = meshes/05_Lucas/elizabeth-canine-lv-endocardium.alg\n")
    file.write("main_function = initialize_grid_with_custom_mesh\n")
    if (use_cluster):
        file.write("library_file = /home/berg/MonoAlg3D_C/shared_libs/libdefault_domains.so.so\n")    
    file.write("\n")

def write_monoalg_purkinje_section (file,seed,rand_offset,algorithm_number,use_cluster=False):
    file.write("[purkinje]\n")
    file.write("name = Elizabeth LV Canine Purkinje\n")
    file.write("start_discretization = 100.0\n")
    file.write("rpmj=10.0\n")
    file.write("pmj_scale=5000.0\n")
    file.write("asymm_ratio=1\n")
    file.write("nmin_pmj=10\n")
    file.write("nmax_pmj=30\n")
    file.write("retro_propagation=true\n")
    file.write("network_file=networks/03_Lucas/01_SRN/03_CO_Activation_Time/seed:%u.vtk\n" % (seed))
    file.write("pmj_location_file=networks/03_Lucas/01_SRN/01_Gold_Standart/elizabeth_pmj_full_um.vtk\n")
    file.write("main_function = initialize_purkinje_coupling_with_custom_mesh\n")
    file.write("library_file = shared_libs/libdefault_purkinje.so\n")
    if (use_cluster):
        file.write("library_file = /home/berg/MonoAlg3D_C/shared_libs/libdefault_purkinje.so\n")
    file.write("\n")

def write_monoalg_purkinje_ode_solver_section (file,use_cluster=False):
    file.write("[purkinje_ode_solver]\n")
    file.write("dt = 0.02\n")
    file.write("use_gpu = no\n")
    file.write("gpu_id = 0\n")
    file.write("library_file = shared_libs/libstewart_aslanidi_noble_2009.so\n")
    if (use_cluster):
        file.write("library_file = /home/berg/MonoAlg3D_C/shared_libs/libstewart_aslanidi_noble_2009.so\n")
    file.write("\n")

def write_monoalg_ode_solver_section (file,use_cluster=False):
    file.write("[ode_solver]\n")
    file.write("dt = 0.02\n")
    file.write("use_gpu = yes\n")
    file.write("gpu_id = 0\n")
    file.write("library_file = shared_libs/libten_tusscher_2006.so\n")
    if (use_cluster):
        file.write("library_file = /home/berg/MonoAlg3D_C/shared_libs/libten_tusscher_2006.so\n")
    file.write("\n")

def write_monoalg_stimulus_section (file,use_cluster=False):
    file.write("[purkinje_stim_his]\n")
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
    #filename = "monoalg3d/co_min:length_seed-%u_nterm:650.ini" % (seed)
    filename = "monoalg3d/co_min:at_seed-%u_nterm:650.ini" % (seed)
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
