#ifndef GRAPH_H
#define GRAPH_H

#include <cstdio>
#include <cstdlib>
#include <cstdbool>
#include <cstdint>
#include <cstring>
#include <vector>
#include <queue>
#include <string>

#include "../utils/utils.h"
#include "../reader/reader.h"

#include "pqueue.h"

// --------------------------------------------------------------------------------
// 'PQUEUE' USER DEFINED STRUCTURES AND FUNCTIONS
typedef struct node_t
{
	pqueue_pri_t pri;
	uint32_t    val;
	size_t pos;
} node_t;


static int
cmp_pri(pqueue_pri_t next, pqueue_pri_t curr)
{
	//return (next < curr);         // equivalent to std::less<int>()
    return (next > curr);           // equivalent to std::greater<int>()
}


static pqueue_pri_t
get_pri(void *a)
{
	return ((node_t *) a)->pri;
}


static void
set_pri(void *a, pqueue_pri_t pri)
{
	((node_t *) a)->pri = pri;
}


static size_t
get_pos(void *a)
{
	return ((node_t *) a)->pos;
}


static void
set_pos(void *a, size_t pos)
{
	((node_t *) a)->pos = pos;
}
// --------------------------------------------------------------------------------

class Edge
{
public:
    uint32_t dest_id;
    double length;
public:
    Edge (const uint32_t dest_id, const double length)
    {
        this->dest_id = dest_id;
        this->length = length;
    }
};

class Node
{
public:
    uint32_t id;
    double x, y, z;
    double lat;
    std::vector<Edge> list_edges;
public:
    Node () { }
    Node (const uint32_t id, const double x,const double y, const double z)
    {
        this->id = id;
        this->x = x;
        this->y = y;
        this->z = z;
    }
    void print ()
    {
        printf("|| [%u] (%g %g %g) {%g} ||",this->id,this->x,this->y,this->z,this->lat);
        for (uint32_t i = 0; i < this->list_edges.size(); i++)
            printf(" --> || [%u] %g || ",this->list_edges[i].dest_id,this->list_edges[i].length);
        printf("\n");
    }
};

class Graph
{
public:
    const double REF_CV = 1900.0;                   // Reference propagation velocity 
    uint32_t total_nodes;                           
    uint32_t total_edges;
    std::vector<Node> list_nodes;                   // List of nodes
    std::vector<double> dist;                       // Shortest distance of each cell from the root point
    std::vector<double> lat;                        // Local activation time of each cell
    std::vector<int> parent;                        // Parent indexes of each cell
    std::vector<uint32_t> terminals_indexes;        // Active PMJ's indexes
public:
    Graph ();
    Graph (std::vector<Point> points, std::vector<Line> lines);
    void depth_first_search (const uint32_t src_id, std::vector<double> &the_segments);
    void dfs (Node u, std::vector<bool> &dfs_visited, std::vector<double> &segments, double &segment_size, uint32_t &flag, uint32_t &total_segments);
    void dijkstra (const uint32_t src_id);
    void dijkstra_2 (const uint32_t src_id);
    uint32_t get_closest_terminal_point (Point p, const bool is_reference);
    void fill_terminal_indexes (std::vector<Point> pmj_points, const bool is_reference);
    void calculate_activation_time (std::vector<Point> pmj_points, const double cv);
    void compute_activation_times (const double cv);
    void compute_errors (Graph *input, std::string pmj_filename);
    void compute_rmse_rrmse_maxerror (std::vector<double> ref_lat, std::vector<double> aprox_lat,\
                                        std::vector<uint32_t> ref_term_ids, std::vector<uint32_t> aprox_term_ids,\
                                        double &rmse, double &rrmse, double &max_error);
    void compute_min_max_lat (std::vector<double> lat, std::vector<uint32_t> term_ids, double &min_value, double &max_value, uint32_t &min_id, uint32_t &max_id);
    void print ();
    void write_terminals (const char filename[]);
};

double adjust_propagation_velocity (const double ref_cv, const double ref_max_lat, const double aprox_max_lat, double &new_diameter);
double calculate_propagation_velocity (const double diameter);
double calculate_proportion (const double ref_diameter, const double ref_max_lat, const double aprox_max_lat);
double calculate_diameter (const double cv);

#endif