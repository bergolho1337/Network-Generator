#ifndef CIRCLE_LIB_H
#define CIRCLE_LIB_H

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
extern "C" void default_cloud_circle (struct cloud_generator_data *generator,\
                                      struct user_data *config);
extern "C" void spaced_cloud_circle (struct cloud_generator_data *generator,\
                                     struct user_data *config);


// Auxiliary functions
void draw_circle_area (const double radius);
void generate_point_inside_circle (double pos[], const double radius);

#endif