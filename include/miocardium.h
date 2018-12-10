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
	double x, y, z;
	double scalar;
	bool taken;
public:
	Point ();
};

class Miocardium
{
public:
	unsigned int num_cloud_points;
	unsigned int num_terminal_points;
	Point *cloud_points;
	Point *terminal_points;
	double max_xyz[3];
	double min_xyz[3];
	
public:
	Miocardium ();

	void read_cloud_points (const std::string miocardium_filename);
	void read_terminal_points (const std::string terminals_filename);	
	void set_limits ();
	void print ();
};

#endif
