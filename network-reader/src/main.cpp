// ********************************************************************************************************************
// 	NETWORK READER
//	Author: Lucas Berg
//
//	Input arguments:
//		- An input file with network written in Legacy VTK file format 
//	Outputs:
//		- Several informations about the structure of teh network 
//
// ********************************************************************************************************************

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <fstream>

#include "graph/graph.h"

using namespace std;

int main (int argc, char *argv[])
{
	if (argc-1 != 1)
	{
		usage(argv[0]);
		exit(EXIT_FAILURE);
	}
	
	Graph *g1 = new Graph(argv[1]);
	//g1->print();

	g1->depth_first_search();

	g1->breadth_first_search();
	
	//g1->write_VTK("outputs/graph-1.vtk");
	
	//g1->remove_edge_graph(4,5);
	//g1->write_VTK("output/graph-1-2.vtk");

	//g2->remove_edge_graph(1,3);
	//g2->write_VTK("output/graph-2-2.vtk");

	//g2->remove_node_graph(1);
	//g2->write_VTK("output/graph-2-3.vtk");

	delete g1;

	return 0;  
}
