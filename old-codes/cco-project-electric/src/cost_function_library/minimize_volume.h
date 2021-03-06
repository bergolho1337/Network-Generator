#ifndef MINIMIZE_VOLUME_LIB_H
#define MINIMIZE_VOLUME_LIB_H

#include <cstdio>
#include <cstdlib>
#include <cmath>

#include "../options/cost_function_config.h"
#include "../options/local_optimization_config.h"
#include "../cco/cco.h"

// Cost functions implementations
extern "C" struct segment_node* minimize_tree_volume_default (struct cco_network *the_network,\
                                        struct cost_function_config *config,\
                                        struct local_optimization_config *local_opt_config,\
                                        const double new_pos[],\
                                        const std::vector<struct segment_node*> feasible_segments);

extern "C" struct segment_node* minimize_tree_volume_with_assymetric_restriction (struct cco_network *the_network,\
                                        struct cost_function_config *config,\
                                        struct local_optimization_config *local_opt_config,\
                                        const double new_pos[],\
                                        const std::vector<struct segment_node*> feasible_segments);

extern "C" struct segment_node* minimize_tree_volume_with_angle_restriction (struct cco_network *the_network,\
                                        struct cost_function_config *config,\
                                        struct local_optimization_config *local_opt_config,\
                                        const double new_pos[],\
                                        const std::vector<struct segment_node*> feasible_segments);

extern "C" struct segment_node* minimize_tree_volume_with_length_restriction (struct cco_network *the_network,\
                                        struct cost_function_config *config,\
                                        struct local_optimization_config *local_opt_config,\
                                        const double new_pos[],\
                                        const std::vector<struct segment_node*> feasible_segments);

extern "C" struct segment_node* minimize_tree_volume_with_angle_and_length_restriction (struct cco_network *the_network,\
                                        struct cost_function_config *config,\
                                        struct local_optimization_config *local_opt_config,\
                                        const double new_pos[],\
                                        const std::vector<struct segment_node*> feasible_segments);

extern "C" struct segment_node* minimize_tree_volume_with_level_penalty (struct cco_network *the_network,\
                                        struct cost_function_config *config,\
                                        struct local_optimization_config *local_opt_config,\
                                        const double new_pos[],\
                                        const std::vector<struct segment_node*> feasible_segments);

extern "C" struct segment_node* minimize_tree_volume_with_level_penalty_and_angle_restriction (struct cco_network *the_network,\
                                        struct cost_function_config *config,\
                                        struct local_optimization_config *local_opt_config,\
                                        const double new_pos[],\
                                        const std::vector<struct segment_node*> feasible_segments);


// Auxiliary functions

#endif