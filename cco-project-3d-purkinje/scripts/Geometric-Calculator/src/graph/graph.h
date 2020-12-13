#ifndef GRAPH_H
#define GRAPH_H

#include <cstdio>
#include <cstdlib>
#include <cstdbool>
#include <cstdint>
#include <cstring>
#include <vector>

#include "../utils/utils.h"

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
    uint32_t total_nodes;
    uint32_t total_edges;
    std::vector<Node> list_nodes;
    std::vector<uint32_t> terminal_indexes;
public:
    Graph ();
    Graph (std::vector<Point> points, std::vector<Line> lines);
    void depth_first_search (const uint32_t src_id, std::vector<double> &the_segments);
    void dfs (Node u, std::vector<bool> &dfs_visited, std::vector<double> &segments, double &segment_size, uint32_t &flag, uint32_t &total_segments);
    void get_edges (std::vector<double> &the_edges);
    void get_bifurcation_angles (std::vector<double> &the_angles);
    void get_segments (std::vector<double> &the_segments);
    uint32_t get_closest_point (Point p);
    void print ();
    void write_info ();
};

#endif