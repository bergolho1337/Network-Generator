#ifndef GENERATOR_H
#define GENERATOR_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <dlfcn.h>

#include <vtkRegularPolygonSource.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkSphereSource.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkDoubleArray.h>
#include <vtkAppendPolyData.h>
#include <vtkCylinderSource.h>
#include <vtkMath.h>
#include <vtkLine.h>

#include "../random/random.h"
#include "../utils/utils.h"
#include "../config/config.h"

// Flag to allow a function export ...
#define EXPORT_FN

// Forward declaration
class Cloud_Generator;

// This works like a function pointer macro ...
#define SET_CLOUD_GENERATOR(name) EXPORT_FN void name(Cloud_Generator *generator,\
                                                      User_Data *config,\
                                                      Random_Generator *random_generator)
typedef SET_CLOUD_GENERATOR(set_cloud_generator_fn);

class Cloud_Generator
{
public:
    uint32_t num_points;                 // Number of points

    std::vector<Point> *points;          // Reference to the cloud points array

    void *handle;                        // Handle to the library that solves a linear system
    set_cloud_generator_fn *function;   // Pointer to the generator function

public:
    Cloud_Generator (User_Data *config);
    ~Cloud_Generator ();
    void set_generator_function (User_Data *config);
    void write_to_pts ();
    void write_to_vtk ();
    void write_to_vtp ();
};

struct cloud_generator_data
{
    uint32_t num_points;                 // Number of points

    std::vector<Point> *points;          // Reference to the cloud points array

    void *handle;                        // Handle to the library that solves a linear system
    set_cloud_generator_fn *function;   // Pointer to the generator function
};

// Constructor and destructor
struct cloud_generator_data* new_cloud_generator (struct user_data *config);
void free_cloud_generator (struct cloud_generator_data *g);

// Auxiliary functions 
void set_generator_function (struct cloud_generator_data *generator, struct user_data *config);
void write_to_pts (struct cloud_generator_data *generator);
void write_to_vtp (struct cloud_generator_data *generator);
void write_to_vtk (struct cloud_generator_data *generator);

#endif
