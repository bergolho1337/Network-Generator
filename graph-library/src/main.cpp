// ********************************************************************************************************************
// 	GRAPH LIBRARY
//	Author: Lucas Berg
//
//	Input arguments:
//		- An input file with the configuration of the graph
//	Outputs:
//		- An output VTK file with the visualization of teh graph
//
// ********************************************************************************************************************

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <fstream>

#include "../include/graph.h"

using namespace std;

int main (int argc, char *argv[])
{
	if (argc-1 != 1)
	{
		usage(argv[0]);
		exit(EXIT_FAILURE);
	}
	
	Graph *g1 = new Graph(argv[1]);
	Graph *g2 = new Graph();

	double pos1[3] = {0.0,0.0,0.0};
	double pos2[3] = {0.0,100.0,0.0};
	double pos3[3] = {-200.0,200.0,0.0};
	double pos4[3] = {200.0,200.0,0.0};
	double pos5[3] = {0.0,400.0,0.0};

	g2->insert_node_graph(pos1);
	g2->insert_node_graph(pos2);
	g2->insert_node_graph(pos3);
	g2->insert_node_graph(pos4);
	g2->insert_node_graph(pos5);

	g2->insert_edge_graph(0,1);
	g2->insert_edge_graph(1,2);
	g2->insert_edge_graph(1,3);
	g2->insert_edge_graph(1,4);

	//g1->print();
	//g2->print();

	g1->write_VTK("output/graph-1.vtk");
	g2->write_VTK("output/graph-2.vtk");

	//g1->remove_edge_graph(4,5);
	//g1->write_VTK("output/graph-1-2.vtk");

	//g2->remove_edge_graph(1,3);
	//g2->write_VTK("output/graph-2-2.vtk");

	//g2->remove_node_graph(1);
	//g2->write_VTK("output/graph-2-3.vtk");

	delete g1;
	delete g2;

	return 0;  
}
