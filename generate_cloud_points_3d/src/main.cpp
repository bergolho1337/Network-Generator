// Author: Lucas Berg

#include "config/config.h"
#include "utils/utils.h"
#include "generator/generator.h"
#include "random/random.h"

using namespace std;


int main (int argc, char *argv[])
{
    if (argc-1 != 1)
    {
        usage(argv[0]);
        exit (EXIT_FAILURE);
    }

    // Configure the random number generator
    const uint32_t seed = 1;
    Random_Generator *random_generator = new Random_Generator(seed);
    random_generator->generate();
    //random_generator->print();

    // Read user input
    User_Data *config = new User_Data();
    config->read(argv[1]);
    config->print();

    Cloud_Generator *generator = new Cloud_Generator(config);

    // Call generator function
    generator->function(generator,config,random_generator);

    // Write the points to file
    generator->write_to_pts();
    generator->write_to_vtk();
    //generator->write_to_vtp();

    // Clean memory
    delete generator;
    delete config;
    delete random_generator;
    
    return 0;
}
