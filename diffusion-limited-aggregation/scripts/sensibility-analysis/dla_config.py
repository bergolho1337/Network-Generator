NUMBER_OF_ITERATIONS = 5000
NUMBER_OF_WALKERS = 3000
ROOT_X = -0.004
ROOT_Y = 0.018
ROOT_Z = 0.0251292

WALKER_RADIUS = 0.0003
WALKER_MESH_FILENAME = "meshes/02_Lucas/elizabeth_exterior_guided_LV.stl"
WALKER_MAP_FILENAME = "maps/02_Lucas/elizabeth_guided_LV_mapping.txt"
WALKER_LIBRARY_NAME = "shared-libs/libcustom_walker.so"

def write_dla_config (seed):
    filename = "dla/dla_seed-%u.ini" % (seed)
    file = open(filename,"w")
    
    # Write [main] section
    file.write("[main]\n")
    file.write("number_of_iterations = %u\n" % NUMBER_OF_ITERATIONS)
    file.write("number_of_walker = %u\n" % NUMBER_OF_WALKERS)
    file.write("root_pos_x = %g\n" % ROOT_X)
    file.write("root_pos_y = %g\n" % ROOT_Y)
    file.write("root_pos_z = %g\n" % ROOT_Z)
    file.write("seed = %u\n" % seed)

    file.write("\n")

    # Write [save_result] section
    file.write("[save_result]\n")
    file.write("output_dir = outputs/01_SRN/03_DLA/seed:%u\n" % (seed))
    file.write("\n")

    # Write [walker] section
    file.write("[walker]\n")
    file.write("use_respawn = false\n")
    file.write("walker_radius = %g\n" % WALKER_RADIUS)
    file.write("mesh_filename = %s\n" % WALKER_MESH_FILENAME)
    file.write("map_filename = %s\n" % WALKER_MAP_FILENAME)
    file.write("library_name = %s\n" % WALKER_LIBRARY_NAME)
    file.write("\n")

    file.close()