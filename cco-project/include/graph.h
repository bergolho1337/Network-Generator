#ifndef GRAPH_H
#define GRAPH_H

#include <cstdio>
#include <cstdlib>
#include <cmath>

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

#endif