//
// Created by bergolho on 03/03/19.
//

#ifndef WALKER_CONFIG_H
#define WALKER_CONFIG_H

#include <cstdio>
#include <cstdlib>
#include <cstdbool>
#include <cstdint>
#include <cstring>
#include <cassert>
#include <dlfcn.h>

#include <map>
#include <string>

#include "../tree/tree.h"

struct walker_config;

// Flag to allow a function export ...
#define EXPORT_FN

#define SET_WALKER_MOVE_FUNCTION(name) EXPORT_FN void name(struct walker *the_walker, struct user_options *the_options)
typedef SET_WALKER_MOVE_FUNCTION(set_walker_move_function_fn);

#define SET_WALKER_RESPAWN_FUNCTION(name) EXPORT_FN void name(double pos[])
typedef SET_WALKER_RESPAWN_FUNCTION(set_walker_respawn_function_fn);

struct walker_config
{
    void *handle;

    char *library_name;
    std::map<std::string,double> *params;    // Parameters of the cost function

    set_walker_move_function_fn *move_function;
    set_walker_respawn_function_fn *respawn_function;
};

struct walker_config* new_walker_config ();
void free_walker_config (struct walker_config *the_walker_config);

void set_walker_functions (struct walker_config *the_walker_config); 

void print_walker_config (struct walker_config *the_walker_config);

#endif