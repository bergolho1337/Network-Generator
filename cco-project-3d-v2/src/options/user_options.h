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
#include "../ini_parser/ini.h"
#include "../ini_parser/ini_file_sections.h"
#include "../single_file_libraries/stb_ds.h"

#include "cost_function_config.h"
#include "local_optimization_config.h"
#include "pruning_config.h"

static const uint32_t MAX_FILENAME_SIZE = 200;

struct user_options
{
    uint32_t seed;
    uint32_t max_rand_offset;
    
    uint32_t n_term;
    double q_perf;
    double p_perf;
    double p_term;
    double v_perf;
    double gamma;

    double root_pos[3];

    bool use_cloud_points;
    char *cloud_points_filename;

    bool use_obstacle;
    char *obstacle_filename;

    bool use_pmj_location;
    char *pmj_location_filename;

    bool use_lat;
    char *lat_filename;

    bool use_only_murray;
    double start_radius;

    struct cost_function_config *config;
    
    bool use_local_optimization;
    struct local_optimization_config *local_opt_config;

    char *output_dir;

    bool use_pruning;
    struct pruning_config *pruning_config;

    bool use_initial_network;
    char *initial_network_filename;

};

struct user_options* new_user_options (int argc, char *argv[]);
void free_user_options (struct user_options *options);

void read_config_file (struct user_options *options, const char filename[]);

int parse_config_file(void *user, const char *section, const char *name, const char *value);

void print_user_options (struct user_options *options);

#endif