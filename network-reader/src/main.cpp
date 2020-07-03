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
#include "pmj/pmj.h"

using namespace std;

void link_pmjs_to_graph(Graph *g, vector<PMJ> pmjs)
{
	for (uint32_t i = 0; i < pmjs.size(); i++)
	{
		double pos1[3];
		pos1[0] = pmjs[i].pos[0];
		pos1[1] = pmjs[i].pos[1];
		pos1[2] = pmjs[i].pos[2];

		Node *min_node = NULL;
		double min_dist = __DBL_MAX__;

		Node *tmp = g->get_list_nodes();
		while (tmp != NULL)
		{
			double pos2[3];
			pos2[0] = tmp->x;
			pos2[1] = tmp->y;
			pos2[2] = tmp->z;
			double dist = sqrt(pow(pos1[0]-pos2[0],2)+pow(pos1[1]-pos2[1],2)+pow(pos1[2]-pos2[2],2));

			if (dist < min_dist)
			{
				min_dist = dist;
				min_node = tmp;
			}

			tmp = tmp->next;
		}
		// Add a new node to the graph as PMJ
		g->insert_node_graph(pos1);
		g->insert_edge_graph(min_node->index,g->get_last_node()->index);
		g->get_last_node()->is_pmj = true;
	}

	g->write_VTK("outputs/pmj_graph.vtk");
}

int main (int argc, char *argv[])
{
	if (argc-1 != 1)
	{
		usage(argv[0]);
		exit(EXIT_FAILURE);
	}
	
	Graph *g1 = new Graph(argv[1]);
	//g1->print();

	//g1->check_duplicates();
	//g1->remove_node_graph(121);

	//g1->depth_first_search();

	//g1->breadth_first_search();

	//g1->write_VTK("outputs/elizabeth_LV.vtk");
	
	//g1->remove_edge_graph(4,5);
	//g1->write_VTK("output/graph-1-2.vtk");

	//g2->remove_edge_graph(1,3);
	//g2->write_VTK("output/graph-2-2.vtk");

	//g2->remove_node_graph(1);
	//g2->write_VTK("output/graph-2-3.vtk");

	//g1->write_network_info();

	// ==============================================================
	// PMJ
	vector<PMJ> pmjs;
	read_pmjs_from_file(pmjs);

	link_pmjs_to_graph(g1,pmjs);

	//g1->write_pmj_config_file("outputs/pmj_config.txt");

	delete g1;

	return 0;  
}
