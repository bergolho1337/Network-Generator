// Author: Lucas Berg

#include "config/config.h"
#include "utils/utils.h"
#include "generator/generator.h"

using namespace std;


int main (int argc, char *argv[])
{
    if (argc-1 != 1)
    {
        usage(argv[0]);
        exit (EXIT_FAILURE);
    }

    // Read user input
    struct user_data *config = new_user_data();

    read_user_input(config,argv[1]);
    print_user_input(config);

    struct cloud_generator_data *generator = new_cloud_generator(config);

    // Call generator function
    generator->function(generator,config);

    // Write the points to file
    //write_to_vtp(generator);
    write_to_vtk(generator);
    write_to_txt(generator);

    free_cloud_generator(generator);
    free_user_data(config);

    return 0;
}
