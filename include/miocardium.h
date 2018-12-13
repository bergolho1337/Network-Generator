#ifndef MIOCARDIUM_H
#define MIOCARDIUM_H

#include <iostream>
#include <fstream>
#include <string>
#include <cstdio>
#include <cstdlib>
#include <cmath>

#include "graph.h"

// Compare function for the C++STL QuickSort
int compare (const void *a, const void *b);

class Point
{
public:
	double x, y, z;				// Coordinates (x,y,z)
	double scalar;				// Value
	bool taken;				// Flag in the case this Point is a terminal
public:
	Point ();
};

class Miocardium
{
public:
	unsigned int num_cloud_points;		// Number of cloud points from the tissue surface
	unsigned int num_terminal_points;	// Number of terminal points from the tissue surface
	Point *cloud_points;			// Reference to the array of cloud points
	Point *terminal_points;			// Reference to the array of terminal points
	double max_xyz[3];			// Minimum limits for the domain (x,y,z)
	double min_xyz[3];			// Maximum limits for the domain (x,y,z)
	
public:
	Miocardium ();

	void read_cloud_points (const std::string miocardium_filename);
	void read_terminal_points (const std::string terminals_filename);	
	void set_limits ();
	void print ();
};

#endif
