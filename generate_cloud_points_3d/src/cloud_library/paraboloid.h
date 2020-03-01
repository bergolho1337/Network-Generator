#ifndef PARABOLOID_LIB_H
#define PARABOLOID_LIB_H

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

#include <vtkQuadric.h>
#include <vtkSampleFunction.h>
#include <vtkContourFilter.h>
#include <vtkOutlineFilter.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkPolyDataWriter.h>
#include <vtkSTLWriter.h>

#include "../generator/generator.h"

// Cost functions implementations
extern "C" void default_paraboloid_cloud (Cloud_Generator *generator,\
                                      User_Data *config,\
                                      Random_Generator *random_generator);


// Auxiliary functions
void draw_default_paraboloid_volume (const double a, const double b, const double side_length);

#endif