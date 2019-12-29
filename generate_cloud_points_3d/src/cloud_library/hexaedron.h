#ifndef HEXAEDRON_LIB_H
#define HEXAEDRON_LIB_H

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
extern "C" void default_hexaedron_cloud (Cloud_Generator *generator,\
                                      User_Data *config,\
                                      Random_Generator *random_generator);


// Auxiliary functions
void draw_default_hexaedron_volume (const double radius);

#endif