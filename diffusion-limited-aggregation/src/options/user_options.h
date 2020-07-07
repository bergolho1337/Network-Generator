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

#include "../ini_parser/ini.h"
#include "../ini_parser/ini_file_sections.h"
#include "../utils/utils.h"

#include "walker_config.h"

struct user_options
{
    uint32_t max_num_iter;
    uint32_t max_num_walkers;
    uint32_t seed;

    double root_pos[3];

    bool use_initial_network;
    char *initial_network_filename;

    struct walker_config *walker_config;
};

struct user_options* new_user_options (int argc, char *argv[]);
void free_user_options (struct user_options *the_options);

void read_config_file (struct user_options *the_options, const char filename[]);
int parse_config_file(void *user, const char *section, const char *name, const char *value);

void print_user_options (struct user_options *the_options);

#endif