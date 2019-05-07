// Author: Lucas Berg

#include <iostream>
#include <cstdio>
#include <cstdlib>

#include "cco/cco.h"
#include "options/user_options.h"
#include "test/test.h"

using namespace std;

int main (int argc, char *argv[])
{   
    if (argc-1 != 1)
    {
        usage(argv[0]);
        exit(EXIT_FAILURE);   
    }

    struct user_options *options = new_user_options(argc,argv);
    struct cco_network *the_network = new_cco_network(options);
    
    grow_tree(the_network,options);
    write_to_vtk(the_network);

    free_cco_network(the_network);
    free_user_options(options);

    return 0;
}