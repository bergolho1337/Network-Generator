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

#define PRINT_LINE "---------------------------------------------------------------------------------"

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
    Edge (const int i, const double weight, Node *destination);

public:
	int index;			// Identifier of the destination node
	double length;		// Size of the edge, euclidean distance
	Edge *next;			// Pointer to the next Edge
	Edge *previous;		// Pointer to the previous Edge
	Node *dest;			// Pointer to the destination Node
};

// =============================================================================================================
// =============================================================================================================
// Strucuture for a Node of the graph
class Node
{
public:
    Node (const int id, const double pos[]);

public:
	int index;			// Identifier of the Node 
	double x, y, z;		// Coordinates (x,y,z)
	int num_edges;		// Number of edges of the Node
	bool is_pmj;		// Flag for a Purkinje-Muscle-Junction node
	Node *next;			// Pointer to the next Node
	Node *previous;		// Pointer to the previous Node

	Edge *list_edges;	// Pointer to the list of Edges
	Edge *last_edge;	// Pointer to the last Edge of the list
};
// =============================================================================================================
// =============================================================================================================
// Structure of the Graph
class Graph
{
public:
    Graph ();
	Graph (const char filename[]);
	~Graph ();

    void print ();
	void write_network_info ();
	double calculate_branch_size (const int parents[], const int ref_index, int &level);
	void check_duplicates ();
	void write_VTK (const char filename[]);
	void write_pmj_config_file (const char filename[]);
	void write_longest_segment (const int parents[], const int ref_index);
	//void printterm ();
    void error (const char msg[]);
	void depth_first_search ();
	void breadth_first_search ();
	//void dijkstra (int s);
	// Inline
	int get_total_nodes () { return total_nodes; }
	int get_total_edges () { return total_edges; }
	Node* get_list_nodes () { return list_nodes; }
	Node* get_last_node () { return last_node; }
	//double* get_dist () { return dist; }
	void build_unitary_vector (Node *u, Node *v, double d[]);
	double calc_angle_between_vectors (const double u[], const double v[]);

	void insert_node_graph (const double pos[]);
	void insert_edge_graph (const int id_1, const int id_2);
	void remove_edge_graph (const int id_1, const int id_2);
	void remove_all_edges (const int id);
	void remove_node_graph (const int id);
	Node* search_node (const int id);

private:
	Node *list_nodes;			// Pointer to the list of Nodes
	Node *last_node;			// Pointer to the last Node of the list
	int total_nodes;			// Total number of Nodes
	int total_edges;			// Total number of Edges
	//double *dist;				// Distance from the source node to all the others

	//bool is_duplicate (const double pos[]);

	void dfs (Node *u, vector<int> &dfs_num);
	void free_list_nodes ();
	void free_list_edges (Node *node);
	void read_graph_from_vtk (const char filename[]);
	void read_graph_from_txt (const char filename[]);

};

// =============================================================================================================
// =============================================================================================================
// Auxiliary functions
void usage (const char pname[]);
void print_stars (const int number);
double calc_norm (double x1, double y1, double z1, double x2, double y2, double z2);
bool check_file_extension (const char filename[], const char extension_name[]);
bool is_terminal (Node *u);
bool is_bifurcation (Node *u);
// =============================================================================================================

#endif
