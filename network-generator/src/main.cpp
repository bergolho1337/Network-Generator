// ********************************************************************************************************************
// 	NETWORK GENERATOR
//	Author: Lucas Berg
//
//	Input arguments:
//		- A configuration file with ".ini" extension (there are sample of these files in the "/examples" folder)
//		- In this configuration file the user must supply:
//			- The cloud of points that represents the tissue where the network will be generated over
//			- The terminal points of the network
//	Outputs:
//		- A Purkinje network that will be stored in the folder "/networks"
//
//---------------------------------------------------------------------------------------------------------------------	
//
// This program builds several types of Purkinje network using different techniques:
// 	L-System
//	CCO
//	Phase Field
// ********************************************************************************************************************

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <omp.h>

#include "../include/generator.h"
#include "../include/timer.h"
#include "../include/config.h"
#include "../include/utils.h"

using namespace std;

int main (int argc, char *argv[])
{
	if (argc-1 != 1)
	{
		Usage(argv[0]);
		exit(EXIT_FAILURE);
	}

	// Allocate all the structures 
	User_Options *options = new User_Options(argc,argv);

	Network_Generator *generator = new Network_Generator(options);

	// Pos-Process the structures here
	// ...

	// Free all the structures
	delete generator;
	delete options;

	
	return 0;  
}
