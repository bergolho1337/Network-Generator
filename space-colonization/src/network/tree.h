#ifndef TREE_H
#define TREE_H

#include <cstdio>
#include <cstdlib>
#include <cstdbool>
#include <cstdint>
#include <cmath>
#include <ctime>
#include <vector>

#include <vtkRegularPolygonSource.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkSphereSource.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkHexahedron.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkDoubleArray.h>
#include <vtkAppendPolyData.h>
#include <vtkCylinderSource.h>
#include <vtkMath.h>
#include <vtkLine.h>

#include "../utils/utils.h"
#include "../options/user_options.h"

#include "leaf.h"
#include "branch.h"
#include "point.h"

const uint32_t WIDTH = 400;
const uint32_t HEIGHT = 400;
const uint32_t MAX_NUMBER_LEAVES = 2000;
const double MIN_DISTANCE = 0.00075;
const double MAX_DISTANCE = 0.001;

class Tree
{
public:
    uint32_t max_number_iterations;
    double min_distance;
    double max_distance;
    double segment_length;

    double root_pos[3];
    double root_dir[3];

    std::vector<Leaf> the_leaves;
    std::vector<Branch> the_branches;

    std::vector<Point> the_cloud_points;
public:
    Tree ();
    Tree (User_Options *options);
    void read_cloud_points (const char filename[]);
    void make_root ();
    void create_leaves (const double leaves_offset);
    void grow_network ();
    void write_branches (uint32_t iter);
    void write_leaves (uint32_t iter);

    void generate();

    void print_cloud_points ();
};

void draw_canvas ();

#endif
