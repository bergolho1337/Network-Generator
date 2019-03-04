#ifndef GRAPH_H_
#define GRAPH_H_

#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <string>
#include <vector>
#include <queue>
#include <algorithm>
#include <fstream>

using namespace std;

const double INF = __DBL_MAX__;				// The value of infinity is set to be the maximum value of a 'double'
const double TOLERANCE_DUPLICATE = 1.0e-20;		// Distance tolerance to consider two Nodes equal

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
	int id;				// Identifier of the destination node
	double w;		    	// Size of the edge, euclidean distance
	Edge *next;			// Pointer to the next Edge
	Node *dest;			// Pointer to the destination Node
};

// =============================================================================================================
// =============================================================================================================
// Strucuture for a Node of the graph
class Node
{
public:
    Node (const int i, const double pos[], const double d[]);

public:
	int id;				// Identifier of the Node 
	double x, y, z;			// Coordinates (x,y,z)
	double d_ori[3];		// Original direction of the growth
	bool is_terminal;		// Flag that identify Node as a terminal or not	
	int num_edges;			// Number of edges of the Node
	Node *next;			// Pointer to the next Node
	Edge *list_edges;		// Pointer to the list of Edges
};
// =============================================================================================================
// =============================================================================================================
// Structure of the Graph
class Graph
{
public:
    Graph ();
	~Graph ();

    	void print ();
	//void printterm ();
    	void error (const char msg[]);
	void dijkstra (int s);
	// Inline
	int get_total_nodes () { return total_nodes; }
	int get_total_edges () { return total_edges; }
	Node* get_list_nodes () { return list_nodes; }
	Node* get_last_node () { return last_node; }
	double* get_dist () { return dist; }

	Node* insert_node_graph (const double pos[], const Node *prev);
	void insert_edge_graph (const int id_1, const int id_2);
	Node* search_node (int id);
private:
	Node *list_nodes;			// Pointer to the lists of Nodes
	Node *last_node;			// Pointer to the last Node of the list
	int total_nodes;			// Total number of Nodes
	int total_edges;			// Total number of Edges
	double *dist;				// Distance from the source node to all the others

	bool is_duplicate (const double pos[]);
	void calc_original_growth_direction (double d_ori[], const Node *prev,\
					    const double x, const double y, const double z);
	void set_term ();

	void free_list_nodes ();
	void free_list_edges (Node *node);

};

// =============================================================================================================
// =============================================================================================================
// Funcoes auxiliares
double calc_norm (double x1, double y1, double z1, double x2, double y2, double z2);
// =============================================================================================================

#endif
