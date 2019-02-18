// Author: Lucas Berg

#include <iostream>
#include <cstdio>
#include <cstdlib>

//#include "point-list/point-list.h"
//#include "segment-list/segment-list.h"
#include "cco/cco.h"
#include "options/user_options.h"

using namespace std;

int main (int argc, char *argv[])
{   
    if (argc-1 != 5)
    {
        printf("Usage:> %s <Qperf> <pperf> <pterm> <rperf> <Nterm>\n",argv[0]);
        exit(EXIT_FAILURE);   
    }

    struct user_options *options = new_user_options(argc,argv);
    struct cco_network *the_network = new_cco_network(options);

    grow_tree(the_network);
    write_to_vtk(the_network);

    return 0;
}