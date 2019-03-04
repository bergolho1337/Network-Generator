#ifndef CO_H
#define CO_H

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <random>
#include <queue>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cfloat>

#include "graph.h"
#include "miocardium.h"
#include "utils.h"

class CO_Config
{
public:
	// CO parameters goes here
	// ...	

	// Common parameters for tree generation
	double root_point[3];
	std::string miocardium_filename;
	std::string terminals_filename;
	
public:
	CO_Config ();
	void print ();
};

class CO_Generator 
{
private:
	Graph *the_purkinje_network;			// Reference to the Purkinje network
	Miocardium *the_miocardium;			// Reference to the miocardium structure
	double root_point[3];				// Root position

public:
	CO_Generator (CO_Config *config);
	void write_network_to_VTK ();
		
};

#endif
