#ifndef PHASEFIELD_H
#define PHASEFIELD_H

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

class PhaseField_Config
{
public:
	// PhaseField parameters goes here
	// ...	

	// Common parameters for tree generation
	double root_point[3];
	std::string miocardium_filename;
	std::string terminals_filename;
	
public:
	PhaseField_Config ();
	void print ();
};

class PhaseField_Generator 
{
private:
	Graph *the_purkinje_network;			// Reference to the Purkinje network
	Miocardium *the_miocardium;			// Reference to the miocardium structure
	double root_point[3];				// Root position

public:
	PhaseField_Generator (PhaseField_Config *config);
	void write_network_to_VTK ();
		
};

#endif
