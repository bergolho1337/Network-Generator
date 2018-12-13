#ifndef LSYSTEM_H
#define LSYSTEM_H

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

#define MAX_SIZE_RAND_ARRAY 1000			// Max size of the random number array
#define SIGMA_RANDOM_DISTRIBUTION 0.2			// 0.1*l_bra gerou bons resultados

#define TOLERANCE_NEAREST 1.0E-08			// Tolerance to find the nearest point over the miocardium

static const double ANGLE = M_PI / 4.0;			// Rotation angle for the growing rule

class Lsystem_Config
{
public:
	double root_point[3];
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
	Graph *the_purkinje_network;			// Reference to the Purkinje network
	Miocardium *the_miocardium;			// Reference to the miocardium structure
	double root_point[3];				// Root position
	double rand_numbers[MAX_SIZE_RAND_ARRAY];	// Array of random numbers
	queue<Node*> growing_nodes;			// Growing queue of Nodes
	int cont_segments;				// Current number of branches from a segment

	// Functions for the L-System tree generation
	void initialize_root_point (const double root_point[]);
	void initialize_random_array (const double l_bra);
	void make_root (const double l_bra);
	void link_to_miocardium ();
	void grow_network (Lsystem_Config *config);
	void calculate_gradient (Node *ptr, double d_gra[], const double cube_size);
	void grow_branch (Node *gnode, const Lsystem_Config *config, const int branch_type);
	double compute_sobel_filter (const Node *p, double sobel[3][3][3], double d_gra[], const double cube_size);
	bool check_mini_cube (const double x, const double y, const double z,\
			     const double width, const double lenght, const double height);
	bool is_inside_cube (const Node *ptr,\
			     const double x, const double y, const double z,\
			     const double width, const double lenght, const double height);
	double calculate_convolution (double d_gra[], double sobel[3][3][3]);
	void calculate_grow_direction (const Node *gnode, const double d_gra[], const Lsystem_Config *config, const double theta);
	void rotate_direction (double d[], const double u[], const double teta);
	void generate_branch (const Node *gnode, const double d[], const Lsystem_Config *config);
	double calculate_size_branch (const double l_bra);
	int search_most_near (const Point *arr, const unsigned int n, const double x, const double y, const double z);
	void join_terminals ();
	void count_number_terminals ();

	bool check_terminals (const Node *gnode, const double x, const double y, const double z, const double tolerance);
	bool check_collision_tree (const Node *gnode, const double x, const double y, const double z, const double tolerance);	
	bool check_collision_miocardium (const Node *gnode, const double x, const double y, const double z, const double tolerance);
	bool check_limits (const Node *gnode, const double x, const double y, const double z);
	

public:
	Lsystem_Generator (Lsystem_Config *config);
	void write_network_to_VTK ();
		
};

#endif
