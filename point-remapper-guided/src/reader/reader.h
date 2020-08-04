#ifndef READER_H_
#define READER_H_

#include <iostream>

#include "graph.h"
#include "cloud.h"

using namespace std;

//#define RADIUS 0.001
#define RADIUS 0.005


class Reader
{
public:
    Graph *the_network;
    Cloud_Point *the_cloud;
public:
    Reader (int argc, char *argv[]);
    void remap_points_using_graph ();
    void check_cloud_points_inside_node (Node *u);
    void write_remapped_points_to_vtk ();
    void write_points_to_pts ();
    bool is_inside (Point_3d point, const double center[]);
};

#endif
