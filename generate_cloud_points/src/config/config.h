#ifndef CONFIG_H
#define GENERATOR_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <dlfcn.h>

struct user_data
{
    double 
};

// Constructor and destructor
struct user_data* new_user_data ();
void free_user_data (struct user_data *u);

// Auxiliary functions 
//char *getMethodName (const int linear_system_id);
//void setMethodFunction (struct linear_system_data *ls, const int linear_system_id);
//void printLinearSystem (struct linear_system_data *ls);
//void printMatrix (const char *name, const double *A, const int n);
//void printVector (const char *name, const double *v, const int n);

#endif