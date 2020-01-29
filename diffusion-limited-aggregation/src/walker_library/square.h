#ifndef MINIMIZE_VOLUME_LIB_H
#define MINIMIZE_VOLUME_LIB_H

#include <cstdio>
#include <cstdlib>
#include <cmath>

#include "../options/walker_config.h"
#include "../tree/tree.h"

// Functions headers
extern "C" void move (struct walker *the_walker, struct user_options *the_options);
extern "C" void respawn (double pos[]);

// Auxiliary functions
bool is_inside (const double width, const double height, const double x, const double y);

#endif