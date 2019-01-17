// Author: Lucas Berg
// This program tries to implementing the Vascular Tree from the paper 
// "Computer-Optimization of Vascular Trees" from Wolfgang Schreiner and Peter Franz Buxbaum 

#include <iostream>

#include "../include/cco.h"

using namespace std;

int main (int argc, char *argv[])
{
    if (argc-1 != 4)
    {
        printf("Usage:> %s <Qperf> <pperf> <rperf> <Nterm>\n",argv[0]);
        exit(EXIT_FAILURE);   
    }

    User_Options *options = read_user_input(argc,argv);
    //print_user_input(options);

    CCO_Network *the_network = new CCO_Network(options);
    
    the_network->grow_tree();
    
    the_network->write_to_vtk();

    delete the_network;

    return 0;
}