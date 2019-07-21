//
// Created by bergolho on 12/02/19.
//

#ifndef USER_OPTIONS_H
#define USER_OPTIONS_H

#include <cstdio>
#include <cstdlib>
#include <cstdbool>
#include <cstdint>
#include <cstring>

#include "../point-list/point-list.h"
#include "../utils/utils.h"

#include "cost_function_config.h"
#include "local_optimization_config.h"

static const uint32_t MAX_FILENAME_SIZE = 200;

struct user_options
{
    int N_term;
    double Q_perf;
    double p_perf;
    double p_term;
    double r_perf;

    double root_pos[3];

    bool use_cloud_points;
    char cloud_points_filename[MAX_FILENAME_SIZE];

    struct cost_function_config *config;
    
    bool use_local_optimization;
    struct local_optimization_config *local_opt_config;

};

struct user_options* new_user_options (int argc, char *argv[]);
void free_user_options (struct user_options *options);

void read_config_file (struct user_options *options, const char filename[]);

void read_main_section (struct user_options *options, FILE *file);
void read_cloud_points_section (struct user_options *options, FILE *file);
void read_cost_function_section (struct user_options *options, FILE *file);
void read_local_optimization_section (struct user_options *options, FILE *file);

void print_user_options (struct user_options *options);

#endif