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

#include "leaf.h"
#include "branch.h"

const uint32_t WIDTH = 400;
const uint32_t HEIGHT = 400;
const uint32_t MAX_NUMBER_LEAVES = 2000;
const double MIN_DISTANCE = 10;
const double MAX_DISTANCE = 100;

class Tree
{
public:
    std::vector<Leaf> the_leaves;
    std::vector<Branch> the_branches;
public:
    Tree ();
    void make_root ();
    void create_leaves ();
    void grow_network ();
    void write_branches (uint32_t iter);
    void write_leaves (uint32_t iter);
};

void draw_canvas ();

#endif
