def write_monoalg_main_section (file,seed,rand_offset):
    file.write("[main]\n")
    file.write("num_threads = 6\n")
    file.write("dt_pde = 0.02\n")
    file.write("simulation_time = 150.0\n")
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
    file.write("print_rate = 10\n")
    #file.write("output_dir = ./outputs/Oxford/min:terminals_seed:%u\n" % (seed))
    #file.write("output_dir = ./outputs/Oxford/min:total_seed:%u\n" % (seed))
    #file.write("output_dir = ./outputs/Oxford/min:terminals-linked_seed:%u\n" % (seed))
    #file.write("output_dir = ./outputs/Oxford/min:total-linked_seed:%u\n" % (seed))
    #file.write("output_dir = ./outputs/Oxford/min:terminals-linked-soft-pruning_seed:%u\n" % (seed))
    #file.write("output_dir = ./outputs/Oxford/min:terminals-linked-moderate-pruning_seed:%u\n" % (seed))
    file.write("output_dir = ./outputs/Oxford/min:terminals-linked-heavy-pruning_seed:%u\n" % (seed))
    file.write("init_function=init_save_purkinje_coupling_with_activation_times\n")
    file.write("end_function=end_save_purkinje_coupling_with_activation_times\n")
    file.write("main_function=save_purkinje_coupling_with_activation_times\n")
    file.write("activation_threshold_tissue = 0.3\n")
    file.write("activation_threshold_purkinje = 0.3\n")
    file.write("apd_threshold_tissue = 0.1\n")
    file.write("apd_threshold_purkinje = 0.1\n")
    file.write("save_activation_time = true\n")
    file.write("save_apd = true\n")
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
    file.write("sigma_x = 0.0000176\n")
    file.write("sigma_y = 0.0000176\n")
    file.write("sigma_z = 0.0000176\n")
    file.write("sigma_purkinje = 0.001\n")
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
    file.write("name = Plain Mesh\n")
    file.write("num_layers = 1\n")
    file.write("start_dx = 100.0\n")
    file.write("start_dy = 100.0\n")
    file.write("start_dz = 100.0\n")
    file.write("side_length=20000\n")
    file.write("main_function = initialize_grid_with_square_mesh\n")
    if (use_cluster):
        file.write("library_file = /home/berg/MonoAlg3D_C/shared_libs/libdefault_domains.so.so\n")    
    file.write("\n")

def write_monoalg_purkinje_section (file,seed,rand_offset,algorithm_number,use_cluster=False):
    file.write("[purkinje]\n")
    file.write("name = Purkinje\n")
    file.write("start_discretization = 100.0\n")
    file.write("rpmj = 1000.0\n")
    file.write("pmj_scale = 500.0\n")
    file.write("nmin_pmj = 10\n")
    file.write("nmax_pmj = 30\n")
    file.write("retro_propagation = true\n")
    #file.write("network_file = networks/03_Lucas/02_Oxford/01_Basics/01_Minimize_Terminals/seed:%u.vtk\n" % (seed))
    #file.write("network_file = networks/03_Lucas/02_Oxford/01_Basics/02_Minimize_Total/seed:%u.vtk\n" % (seed))
    #file.write("network_file = networks/03_Lucas/02_Oxford/01_Basics/03_Minimize_Terminals_Linked/seed:%u.vtk\n" % (seed))
    #file.write("network_file = networks/03_Lucas/02_Oxford/01_Basics/04_Minimize_Total_Linked/seed:%u.vtk\n" % (seed))
    #file.write("network_file = networks/03_Lucas/02_Oxford/01_Basics/05_Minimize_Terminals_Linked_Soft_Pruning/seed:%u.vtk\n" % (seed))
    #file.write("network_file = networks/03_Lucas/02_Oxford/01_Basics/06_Minimize_Terminals_Linked_Moderate_Pruning/seed:%u.vtk\n" % (seed))
    file.write("network_file = networks/03_Lucas/02_Oxford/01_Basics/07_Minimize_Terminals_Linked_Heavy_Pruning/seed:%u.vtk\n" % (seed))
    file.write("pmj_location_file = networks/03_Lucas/02_Oxford/01_Basics/00_Reference/reference_pmjs.vtk\n")
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
    file.write("library_file = shared_libs/libfhn_mod.so\n")
    if (use_cluster):
        file.write("library_file = /home/berg/MonoAlg3D_C/shared_libs/libfhn_mod.so\n")
    file.write("\n")

def write_monoalg_ode_solver_section (file,use_cluster=False):
    file.write("[ode_solver]\n")
    file.write("dt = 0.02\n")
    file.write("use_gpu = yes\n")
    file.write("gpu_id = 0\n")
    file.write("library_file = shared_libs/libmitchell_shaeffer_2003.so\n")
    if (use_cluster):
        file.write("library_file = /home/berg/MonoAlg3D_C/shared_libs/libmitchell_shaeffer_2003.so\n")
    file.write("\n")

def write_monoalg_stimulus_section (file,use_cluster=False):
    file.write("[purkinje_stim_his]\n")
    file.write("start = 0.0\n")
    file.write("duration = 2.0\n")
    file.write("current = 1.0\n")
    file.write("x_limit = 500\n")
    file.write("main_function = stim_if_x_less_than\n")
    if (use_cluster):
        file.write("library_file = /home/berg/MonoAlg3D_C/shared_libs/libdefault_stimuli.so")
    file.write("\n")


def write_monoalg_config_file (seed,rand_offset,algorithm_number):
    filename = "monoalg3d/min:terminals-linked-heavy_seed-%u.ini" % (seed)
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
