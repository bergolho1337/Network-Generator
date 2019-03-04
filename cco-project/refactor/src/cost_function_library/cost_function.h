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


// Auxiliary functions

#endif