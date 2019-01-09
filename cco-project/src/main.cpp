#include <iostream>

#include "../include/cco.h"

using namespace std;

int main (int argc, char *argv[])
{
    if (argc-1 != 6)
    {
        printf("Usage:> %s <x0> <y0> <Qperf> <pperf> <rperf> <Nterm>\n",argv[0]);
        exit(EXIT_FAILURE);   
    }

    User_Options *options = read_user_input(argc,argv);
    //print_user_input(options);

    Graph *the_network = initialize_graph();

    grow_cco_tree(the_network,options);
    write_graph_to_VTK(the_network);
    //print_graph(the_network);

    //free_network(the_network);

    return 0;
}