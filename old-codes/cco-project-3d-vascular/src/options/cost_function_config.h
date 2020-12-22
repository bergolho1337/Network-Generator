//
// Created by bergolho on 03/03/19.
//

#ifndef COST_FUNCTION_CONFIG_H
#define COST_FUNCTION_CONFIG_H

#include <cstdio>
#include <cstdlib>
#include <cstdbool>
#include <cstdint>
#include <cstring>
#include <dlfcn.h>

#include <map>
#include <string>

#include "../cco/cco.h"

#include "user_options.h"

struct cost_function_config;

// Flag to allow a function export ...
#define EXPORT_FN

// Template for a cost function
#define SET_COST_FUNCTION(name) EXPORT_FN struct segment_node* name(struct cco_network *the_network,\
                                                struct cost_function_config *config,\
                                                struct local_optimization_config *local_opt_config,\
                                                const double new_pos[],\
                                                const std::vector<struct segment_node*> feasible_segments,\
                                                const std::vector<struct face> obstacle_faces,\
                                                const std::vector<struct pmj_lat> lat_points)
typedef SET_COST_FUNCTION(set_cost_function_fn);

struct cost_function_config
{
    void *handle;

    char *function_name;
    char *library_name;
    std::map<std::string,double> *params;    // Parameters of the cost function

    set_cost_function_fn *function;       // Reference to the cost function

};

struct cost_function_config* new_cost_function_config ();
void free_cost_function_config (struct cost_function_config *config);

void set_cost_function (struct cost_function_config *config);

bool get_parameter_value_from_map (std::map<std::string,double> *params,\
                                const std::string key, double *value);

void print_cost_function_config (struct cost_function_config *config);

#endif