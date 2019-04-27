#ifndef LOCAL_OPTIMIZATION_LIB_H
#define LOCAL_OPTIMIZATION_LIB_H

#include <cstdio>
#include <cstdlib>
#include <cmath>

#include "../options/local_optimization_config.h"
#include "../cco/cco.h"

// Local optimization functions implementations
extern "C" void default_local_optimization (struct segment_node *iconn,\
                     struct segment_node *ibiff,\
                     struct segment_node *inew,\
                     std::vector<struct point*> &test_positions);

// Auxiliary functions


#endif