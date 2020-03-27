#ifndef SQUARE_LIB_H
#define SQUARE_LIB_H

#include <cstdio>
#include <cstdlib>
#include <cmath>

#include <vtkRegularPolygonSource.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkHexahedron.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkDoubleArray.h>
#include <vtkMath.h>
#include <vtkLine.h>

#include "../generator/generator.h"

// Cost functions implementations
extern "C" void default_square_cloud (Cloud_Generator *generator,\
                                      User_Data *config,\
                                      Random_Generator *random_generator);


// Auxiliary functions
void draw_default_square_volume (const double side_length);

#endif