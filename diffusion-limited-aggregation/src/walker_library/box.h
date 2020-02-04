#ifndef BOX_LIB_H
#define BOX_LIB_H

#include <cstdio>
#include <cstdlib>
#include <cmath>

#include "../options/walker_config.h"
#include "../tree/tree.h"

// Functions headers
extern "C" void move (struct walker *the_walker, struct user_options *the_options);
extern "C" void respawn (struct walker_config *the_walker_config, double pos[]);
extern "C" void draw (struct walker_config *the_walker_config, const double root_pos[]);

// Auxiliary functions
bool is_inside (const double width, const double height, const double depth,\
                const double x, const double y, const double z);

#endif