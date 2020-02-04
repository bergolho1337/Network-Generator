//
// Created by bergolho on 26/01/20.
//

#ifndef TREE_H
#define TREE_H

#include <cstdio>
#include <cstdlib>
#include <cstdbool>
#include <cstdint>
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

#include "walker.h"
#include "segment.h"

#include "../options/user_options.h"
#include "../utils/stop_watch.h"

struct dla_tree
{
    struct walker_list *point_list;
    struct segment_list *segment_list;
};

struct dla_tree* new_dla_tree ();
void free_dla_tree (struct dla_tree *the_tree);

void grow_tree (struct dla_tree *the_tree, struct user_options *the_options);
void write_to_vtk (struct dla_tree *the_tree);

void print_dla_tree (struct dla_tree *the_tree);

#endif