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
#include "network.h"

#include "../options/user_options.h"
#include "../utils/stop_watch.h"

const uint32_t MAX_NUMBER_OF_NODES = 1000;

struct dla_tree
{
    struct walker_list *point_list;
    struct segment_list *segment_list;
};

struct dla_tree* new_dla_tree ();
void free_dla_tree (struct dla_tree *the_tree);

void make_root_default (struct dla_tree *the_tree, struct user_options *the_options);
void make_root_using_initial_network (struct dla_tree *the_tree, struct user_options *the_options);
void grow_tree (struct dla_tree *the_tree, struct user_options *the_options);

void write_tree_to_vtk (struct dla_tree *the_tree, const char output_dir[]);
void write_current_tree_to_vtk (struct dla_tree *the_tree, const char output_dir[]);
void write_root (struct walker_list *l, const char output_dir[]);

void print_dla_tree (struct dla_tree *the_tree);

#endif