#ifndef GRAPH_H
#define GRAPH_H

#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cassert>
#include <string>
#include <vector>
#include <queue>
#include <algorithm>
#include <fstream>

using namespace std;

// Constant and macros
#define PRINT_LINE "---------------------------------------------------------------------------------"

const double INF = __DBL_MAX__;				// The value of infinity is set to be the maximum value of a 'double'
const double TOLERANCE_DUPLICATE = 1.0e-20;		// Distance tolerance to consider two Nodes equal
const double GAMMA = 2.55;       // Bifurcation expoent
const double ETA = 3.6e-03;      // Viscosity of blood


class Edge;
class Node;

// =============================================================================================================
// =============================================================================================================
// Structure for an Edge of the graph (a.k.a Segment)
class Edge
{
public:
    Edge (const int i, const double weight, Node *destination);

public:
	int index;			// Identifier of the destination node
	double length;		// Size of the edge, euclidean distance
    double radius;      // Radius of the Segment
    double pressure;    // Pressure of the Segment
    double biff_ratio;  // Bifurcation ratio of the Segment

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
	void write_VTK (const char filename[]);
	//void printterm ();
    void error (const char msg[]);
	//void dijkstra (int s);
	// Inline
	int get_total_nodes () { return total_nodes; }
	int get_total_edges () { return total_edges; }
	Node* get_list_nodes () { return list_nodes; }
	Node* get_last_node () { return last_node; }
	//double* get_dist () { return dist; }

	void insert_node_graph (const double pos[]);
	void insert_edge_graph (const int id_1, const int id_2, const double Q_term, const double p_term);
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

	void free_list_nodes ();
	void free_list_edges (Node *node);

};

// =============================================================================================================
// =============================================================================================================
// Auxiliary functions
double calc_norm (double x1, double y1, double z1, double x2, double y2, double z2);
double poiseuille (const double Q_term, const double p_term, const double length);
// =============================================================================================================


/*
struct Edge // Segment
{
    int index;                  // Destination index
    int type;                   // Type of the Segment (0 = parent, 1 = left, 2 = right)
    double length;              // Length of the Segment
    double radius;              // Radius of the Segment
    double pressure;            // Pressure of the Segment
    double biff_ratio;          // Bifurcation ratio of the Segment
    double resistance;          // Relative resistance of the Segment
    double old_length;          // Old length of the segment

    Edge *next;
}typedef Edge;

struct Node 
{
    int index;                  // Source index
    double x, y;                // Coordinate of the point
    Node *parent;               // Pointer to its parent
    Node *left_offspring;       // Pointer to the left offsping
    Node *right_offpring;       // Pointer to the right offspring
    // Node *offsprings[2];     // This should be better ...
    int num_segments;

    Node *next;                 // Next node from the list of nodes
    Edge *segment_list;         // List of edges


}typedef Node;

struct Graph
{
    Node *list_nodes;           // List of nodes
    Node *last_node;            // Pointer to the last node
    int num_nodes;              // Number of nodes
    int num_segments;           // Number of segments
}typedef Graph;

struct Terminal
{
    int node_index;             // Node index of the terminal 
    double pressure_term;       // Pressure of the terminal
    double flow_term;           // Flow of the terminal
}typedef Terminal;

Graph* initialize_graph ();
void insert_node_graph (Graph *g, const double x, const double y);
void insert_edge_graph (Graph *g, const int source, const int destination);
Node* search_node (Graph *g, const int id);
double calc_segment_length (const Node *ptr1, const Node *ptr2);
void print_graph (Graph *g);
void write_graph_to_VTK (Graph *g);
*/

#endif