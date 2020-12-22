#ifndef PRUNING_CONFIG_H
#define PRUNING_CONFIG_H

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <map>
#include <dlfcn.h>
#include <cassert>

struct pruning_config;

// Flag to allow a function export ...
#define EXPORT_FN

// Template for a local optimization function
#define SET_PRUNING_FUNCTION(name) EXPORT_FN double name(struct pruning_config *config, const double level, const double length)
typedef SET_PRUNING_FUNCTION(set_pruning_function_fn);

struct pruning_config
{
    void *handle;

    char *function_name;
    char *library_name;
    std::map<std::string,double> *params;       // Parameters of the pruning function

    set_pruning_function_fn *function;          // Reference to the pruning function

};

struct pruning_config* new_pruning_config ();
void free_pruning_config (struct pruning_config *config);

void set_pruning_function (struct pruning_config *config);

void print_pruning_config (struct pruning_config *config);

// Auxiliary functions

#endif