#ifndef GRAPH_H_
#define GRAPH_H_

#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cassert>
#include <cstring>
#include <string>
#include <vector>
#include <queue>
#include <map>
#include <algorithm>
#include <fstream>

#include "../utils/utils.h"

using namespace std;

const double INF = __DBL_MAX__;				// The value of infinity is set to be the maximum value of a 'double'
const double TOLERANCE_DUPLICATE = 1.0e-20;		// Distance tolerance to consider two Nodes equal
const int DFS_WHITE = -1;
const int DFS_BLACK = 1;

class Edge;
class Node;

// =============================================================================================================
// =============================================================================================================
// Structure for an Edge of the graph
class Edge
{
public:
	Edge () { }
    Edge (const uint32_t dest_id, const double length)
	{
		this->dest_index = dest_id;
		this->length = length;
	}
public:
	uint32_t dest_index;			// Identifier of the destination node
	double length;					// Size of the edge, euclidean distance
};

// =============================================================================================================
// =============================================================================================================
// Strucuture for a Node of the graph
class Node
{
public:
	Node () { }
    Node (const int id, const double pos[])
	{
		this->index = id;
		this->x = pos[0];
		this->y = pos[1];
		this->z = pos[2];
		this->is_pmj = false;
	}
	void set_node (const uint32_t id, const double pos[])
	{
		this->index = id;
		this->x = pos[0];
		this->y = pos[1];
		this->z = pos[2];
		this->is_pmj = false;
	}
public:
	uint32_t index;						// Identifier of the Node 
	double x, y, z;						// Coordinates (x,y,z)
	bool is_pmj;						// Flag for a Purkinje-Muscle-Junction node
	std::vector<Edge> list_edges;		// List of edges
};
// =============================================================================================================
// =============================================================================================================
// Structure of the Graph
class Graph
{
public:
    Graph ();
	Graph (const char filename[]);
	void read_graph_from_vtk (const char filename[]);
	void read_graph_from_txt (const char filename[]);
	void build_unitary_vector (double d[], const uint32_t src_id, const uint32_t dest_id);
	void set_terminals ();

	void depth_first_search (const uint32_t src_id);
	void dfs (const uint32_t src_id, std::vector<bool> &dfs_num);

    void print ();
	void write_terminals ();
	void write_network_info ();
	bool check_duplicates ();
	bool is_bifurcation (Node u);
	
	//void write_network_info ();
	//double calculate_branch_size (const int parents[], const int ref_index, int &level);
	//void check_duplicates ();
	//void write_VTK (const char filename[]);
	//void write_pmj_config_file (const char filename[]);
	//void write_longest_segment (const int parents[], const int ref_index);
	//void print_terminals ();
    //void error (const char msg[]);
	//void depth_first_search ();
	//void breadth_first_search ();
	//void dijkstra (int s);
	//double* get_dist () { return dist; }
	//void build_unitary_vector (Node *u, Node *v, double d[]);
	//double calc_angle_between_vectors (const double u[], const double v[]);

	//void insert_node_graph (const double pos[]);
	//void insert_edge_graph (const int id_1, const int id_2);
	//void remove_edge_graph (const int id_1, const int id_2);
	//void remove_all_edges (const int id);
	//void remove_node_graph (const int id);
	//Node* search_node (const int id);

public:
	std::vector<Node> list_nodes;
	uint32_t number_of_terminals;
	uint32_t total_nodes;
	uint32_t total_edges;

	//double *dist;				// Distance from the source node to all the others

	//bool is_duplicate (const double pos[]);

	//void dfs (Node *u, vector<int> &dfs_num);
	//void free_list_nodes ();
	//void free_list_edges (Node *node);
	//void read_graph_from_vtk (const char filename[]);
	//void read_graph_from_txt (const char filename[]);

};

// =============================================================================================================
// =============================================================================================================
// Auxiliary functions
bool is_terminal (Node *u);
bool is_bifurcation (Node *u);
// =============================================================================================================

#endif
