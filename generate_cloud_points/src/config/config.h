#ifndef CONFIG_H
#define CONFIG_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <dlfcn.h>

#include <string>
#include <map>
#include <fstream>

#include "../utils/utils.h"

struct user_data
{
    u_int32_t num_points;
    double area;
    std::string *function_name;
    std::string *library_name;
    
    std::map<std::string,double> *param;
};

// Constructor and destructor
struct user_data* new_user_data ();
void free_user_data (struct user_data *config);

void read_user_input (struct user_data *config, const char filename[]);
void print_user_input (struct user_data *config);

// Auxiliary functions 
//char *getMethodName (const int linear_system_id);
//void setMethodFunction (struct linear_system_data *ls, const int linear_system_id);
//void printLinearSystem (struct linear_system_data *ls);
//void printMatrix (const char *name, const double *A, const int n);
//void printVector (const char *name, const double *v, const int n);

#endif