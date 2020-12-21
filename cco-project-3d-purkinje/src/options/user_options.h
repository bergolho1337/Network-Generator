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
#include <string>

//#include "../point-list/point-list.h"
//#include "../utils/utils.h"
#include "../ini_parser/ini.h"
#include "../ini_parser/ini_file_sections.h"
#include "../single_file_libraries/stb_ds.h"

#include "cost_function_config.h"
#include "local_optimization_config.h"
#include "pmj_config.h"
//#include "pruning_config.h"

#define PRINT_LINE "====================================================================================="
#define PRINT_DOTS "....................................................................................."

class User_Options
{
public:
    uint32_t n_term;
    uint32_t seed;
    uint32_t max_rand_offset;
    
    double root_pos[3];

    bool use_only_murray;
    double start_radius;
    double gamma;
    double lat_offset;

    bool use_cloud_points;
    std::string cloud_points_filename;

    bool use_obstacle;
    std::string obstacle_filename;

    bool use_pmj_location;
    PMJConfig *pmj_config;

    bool use_initial_network;
    std::string initial_network_filename;

    std::string output_dir;

    CostFunctionConfig *cost_function_config;
    
    bool use_local_optimization;
    LocalOptimizationConfig *local_opt_config;

    bool use_pruning;
    //struct pruning_config *pruning_config;

public:
    User_Options (const char filename[]);
    ~User_Options ();
    void read_config_file (const char filename[]);
    void print ();
};

int parse_config_file(void *user, const char *section, const char *name, const char *value);

#endif