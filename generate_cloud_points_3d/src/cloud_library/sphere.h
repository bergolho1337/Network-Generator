#ifndef SPHERE_LIB_H
#define SPHERE_LIB_H

#include <cstdio>
#include <cstdlib>
#include <cmath>

#include <vtkRegularPolygonSource.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkSphereSource.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkDoubleArray.h>
#include <vtkCylinderSource.h>
#include <vtkMath.h>
#include <vtkLine.h>

#include "../generator/generator.h"

// Cost functions implementations
extern "C" void default_sphere_cloud (Cloud_Generator *generator,\
                                      User_Data *config,\
                                      Random_Generator *random_generator);


// Auxiliary functions
void draw_default_sphere_volume (const double radius);

#endif