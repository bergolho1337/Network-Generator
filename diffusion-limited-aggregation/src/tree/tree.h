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

#include "vtkPoints.h"
#include "vtkPointData.h"
#include "vtkXMLPolyDataReader.h"
#include "vtkPolyData.h"
#include <vtkVersion.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkSmartPointer.h>
#include <vtkDoubleArray.h>
#include <vtkAppendPolyData.h>
#include <vtkSphereSource.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkHexahedron.h>
#include <vtkUnstructuredGrid.h>

#include "walker.h"
#include "segment.h"

const uint32_t MAX_NUMBER_OF_WALKERS = 50;
const uint32_t MAX_NUMBER_OF_ITERATIONS = 1500000;

class Tree
{
public:
    std::vector<Walker*> the_tree;
    std::vector<Segment*> the_segments;
public:
    Tree ();
    void grow ();
    void write_to_vtk (const uint32_t iter);
    void write_to_vtp ();
};

#endif