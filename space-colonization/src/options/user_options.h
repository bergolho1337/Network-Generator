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
    uint32_t max_iterations;

    double initial_length;
    double initial_angle;
    double initial_diameter;

    double angle_decrease_ratio;
    double length_decrease_ratio;
    double diameter_decrease_ratio;

    double root_pos[3];

public:
    User_Options (int argc, char *argv[]);
    
    void read_config_file (const char filename[]);

    void print ();
};

int parse_config_file(void *user, const char *section, const char *name, const char *value);


#endif