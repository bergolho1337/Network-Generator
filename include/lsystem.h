#ifndef LSYSTEM_H
#define LSYSTEM_H

#include <iostream>
#include <fstream>
#include <string>
#include <random>
#include <cstdio>
#include <cstdlib>
#include <cmath>

#include "graph.h"
#include "miocardium.h"

#define MAX_SIZE_RAND_ARRAY 1000			// Max size of the random number array
#define SIGMA_RANDOM_DISTRIBUTION 0.2			// 0.1*l_bra gerou bons resultados

class Lsystem_Config
{
public:
	unsigned int max_grow_iterations;
	double branch_length;
	double w_1;
	double sobel_cube_size;
	double tolerance_tree_collision;
	double tolerance_miocardium_collision;
	double tolerance_terminal_collision;
	std::string miocardium_filename;
	std::string terminals_filename;
	
public:
	Lsystem_Config ();
	void print ();
};

class Lsystem_Generator 
{
private:
	Graph *the_purkinje_network;
	Miocardium *the_miocardium;
	double root_point[3];
	double rand_numbers[MAX_SIZE_RAND_ARRAY];

	void initialize_root_point ();
	void initialize_random_array (const double l_bra);

public:
	Lsystem_Generator (Lsystem_Config *config);
		
};

#endif
