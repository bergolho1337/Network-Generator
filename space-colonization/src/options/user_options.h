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

#include "../utils/utils.h"
#include "../ini_parser/ini.h"
#include "../ini_parser/ini_file_sections.h"

static const uint32_t MAX_FILENAME_SIZE = 200;

class User_Options
{
public:
    uint32_t number_of_iterations;

    double segment_length;
    double root_pos[3];
    double root_dir_pos[3];

    uint32_t leaves_offset;

    double min_distance;
    double max_distance;

    char *cloud_points_filename;

public:
    User_Options (int argc, char *argv[]);
    
    void read_config_file (const char filename[]);

    void print ();
};

int parse_config_file(void *user, const char *section, const char *name, const char *value);


#endif