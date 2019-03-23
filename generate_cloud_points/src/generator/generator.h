#ifndef GENERATOR_H
#define GENERATOR_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <dlfcn.h>

#include "../utils/utils.h"

// Flag to allow a function export ...
#define EXPORT_FN

// This works like a function pointer macro ...
#define SET_CLOUD_GENERATOR(name) EXPORT_FN void name(const double a, const double b)
typedef SET_CLOUD_GENERATOR(set_cloud_generator_fn);

struct cloud_generator_data
{
    void *handle;                        // Handle to the library that solves a linear system
    char *method_name;                   // Name of the generator method
    
    
    set_cloud_generator_fn *generator;   // Pointer to the generator function
};

// Constructor and destructor
struct cloud_generator_data* new_cloud_generator ();
void free_cloud_generator (struct cloud_generator_data *g);

// Auxiliary functions 
//char *getMethodName (const int linear_system_id);
//void setMethodFunction (struct linear_system_data *ls, const int linear_system_id);
//void printLinearSystem (struct linear_system_data *ls);
//void printMatrix (const char *name, const double *A, const int n);
//void printVector (const char *name, const double *v, const int n);

#endif