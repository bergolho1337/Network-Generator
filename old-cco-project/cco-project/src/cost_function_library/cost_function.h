#ifndef COST_FUNCTION_LIB_H
#define COST_FUNCTION_LIB_H

#include <cstdio>
#include <cstdlib>
#include <cmath>

#include "../options/cost_function_config.h"
#include "../cco/cco.h"

// Cost functions implementations
extern "C" struct segment_node* funcao1 (struct cco_network *the_network,\
                    struct cost_function_config *config,\
                    const double new_pos[],\
                    const std::vector<struct segment_node*> feasible_segments);

extern "C" struct segment_node* funcao2 (struct cco_network *the_network,\
                    struct cost_function_config *config,\
                    const double new_pos[],\
                    const std::vector<struct segment_node*> feasible_segments);

extern "C" struct segment_node* funcao3 (struct cco_network *the_network,\
                    struct cost_function_config *config,\
                    const double new_pos[],\
                    const std::vector<struct segment_node*> feasible_segments);

extern "C" struct segment_node* funcao4 (struct cco_network *the_network,\
                    struct cost_function_config *config,\
                    const double new_pos[],\
                    const std::vector<struct segment_node*> feasible_segments);

extern "C" struct segment_node* closest_segment (struct cco_network *the_network,\
                    struct cost_function_config *config,\
                    const double new_pos[],\
                    const std::vector<struct segment_node*> feasible_segments);

extern "C" struct segment_node* closest_segment_with_limit_size (struct cco_network *the_network,\
                    struct cost_function_config *config,\
                    const double new_pos[],\
                    const std::vector<struct segment_node*> feasible_segments);
                
extern "C" struct segment_node* closest_segment_with_angle_restriction (struct cco_network *the_network,\
                    struct cost_function_config *config,\
                    const double new_pos[],\
                    const std::vector<struct segment_node*> feasible_segments);

extern "C" struct segment_node* minimize_tree_volume_default (struct cco_network *the_network,\
                    struct cost_function_config *config,\
                    const double new_pos[],\
                    const std::vector<struct segment_node*> feasible_segments);

extern "C" struct segment_node* minimize_tree_volume_with_assymetric_restriction (struct cco_network *the_network,\
                    struct cost_function_config *config,\
                    const double new_pos[],\
                    const std::vector<struct segment_node*> feasible_segments);

extern "C" struct segment_node* minimize_tree_volume_with_angle_restriction (struct cco_network *the_network,\
                    struct cost_function_config *config,\
                    const double new_pos[],\
                    const std::vector<struct segment_node*> feasible_segments);

extern "C" struct segment_node* maximize_tree_volume (struct cco_network *the_network,\
                    struct cost_function_config *config,\
                    const double new_pos[],\
                    const std::vector<struct segment_node*> feasible_segments);

extern "C" struct segment_node* minimize_tree_volume_with_level_penalty (struct cco_network *the_network,\
                    struct cost_function_config *config,\
                    const double new_pos[],\
                    const std::vector<struct segment_node*> feasible_segments);

extern "C" struct segment_node* minimize_tree_activation_time (struct cco_network *the_network,\
                    struct cost_function_config *config,\
                    const double new_pos[],\
                    const std::vector<struct segment_node*> feasible_segments);

extern "C" struct segment_node* minimize_tree_activation_time_with_angle_restriction (struct cco_network *the_network,\
                    struct cost_function_config *config,\
                    const double new_pos[],\
                    const std::vector<struct segment_node*> feasible_segments);

extern "C" struct segment_node* minimize_tree_activation_time_with_angle_restriction_and_level_restriction (struct cco_network *the_network,\
                    struct cost_function_config *config,\
                    const double new_pos[],\
                    const std::vector<struct segment_node*> feasible_segments);

extern "C" struct segment_node* minimize_custom_function (struct cco_network *the_network,\
                    struct cost_function_config *config,\
                    const double new_pos[],\
                    const std::vector<struct segment_node*> feasible_segments);

extern "C" struct segment_node* minimize_custom_function_with_level_penalty (struct cco_network *the_network,\
                    struct cost_function_config *config,\
                    const double new_pos[],\
                    const std::vector<struct segment_node*> feasible_segments);

// Auxiliary functions

#endif