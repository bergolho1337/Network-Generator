#ifndef READER_H_
#define READER_H_

#include <iostream>

#include "graph.h"
#include "cloud.h"

using namespace std;

// Elizabeth RV config
//#define RADIUS 0.007
//#define SCALE_RATIO 0.25

// Rafa Sebastina config
//#define RADIUS 4.0
//#define SCALE_RATIO 0.001

#define RADIUS 0.01
#define SCALE_RATIO 0.25

class Reader
{
public:
    Graph *the_network;
    Cloud_Point *the_cloud;
public:
    Reader (int argc, char *argv[]);
    void remap_points_using_graph ();
    void depth_first_search (const uint32_t source_index);
    void dfs (Node *u, vector<int> &dfs_num);
    void check_cloud_points_inside_node (Node *u);
    void write_remapped_points_to_vtk ();
    void write_points_to_pts ();
    bool is_inside (Point_3d point, const double center[]);
};

#endif
