#ifndef CLOSEST_SEGMENT_LIB_H
#define CLOSEST_SEGMENT_LIB_H

#include <cstdio>
#include <cstdlib>
#include <cmath>

#include "../options/cost_function_config.h"
#include "../options/local_optimization_config.h"
#include "../cco/cco.h"

// Cost functions implementations
extern "C" struct segment_node* closest_segment (struct cco_network *the_network,\
                                        struct cost_function_config *config,\
                                        struct local_optimization_config *local_opt_config,\
                                        const double new_pos[],\
                                        const std::vector<struct segment_node*> feasible_segments);

extern "C" struct segment_node* closest_segment (struct cco_network *the_network,\
                                        struct cost_function_config *config,\
                                        struct local_optimization_config *local_opt_config,\
                                        const double new_pos[],\
                                        const std::vector<struct segment_node*> feasible_segments);

extern "C" struct segment_node* closest_segment_with_limit_size (struct cco_network *the_network,\
                                        struct cost_function_config *config,\
                                        struct local_optimization_config *local_opt_config,\
                                        const double new_pos[],\
                                        const std::vector<struct segment_node*> feasible_segments);

extern "C" struct segment_node* closest_segment_with_angle_restriction (struct cco_network *the_network,\
                                        struct cost_function_config *config,\
                                        struct local_optimization_config *local_opt_config,\
                                        const double new_pos[],\
                                        const std::vector<struct segment_node*> feasible_segments);

// Auxiliary functions

#endif