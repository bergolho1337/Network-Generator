#ifndef MINIMIZE_ACTIVATION_TIME_LIB_H
#define MINIMIZE_ACTIVATION_TIME_LIB_H

#include <cstdio>
#include <cstdlib>
#include <cmath>

#include "../options/cost_function_config.h"
#include "../options/local_optimization_config.h"
#include "../cco/cco.h"

// Cost functions implementations
extern "C" struct segment_node* minimize_activation_time (struct cco_network *the_network,\
                                        struct cost_function_config *config,\
                                        struct local_optimization_config *local_opt_config,\
                                        const double new_pos[],\
                                        const std::vector<struct segment_node*> feasible_segments,\
                                        const std::vector<struct face> obstacle_faces,\
                                        const std::vector<struct pmj_lat> lat_points);
extern "C" struct segment_node* minimize_activation_time_with_angle_restriction (struct cco_network *the_network,\
                                        struct cost_function_config *config,\
                                        struct local_optimization_config *local_opt_config,\
                                        const double new_pos[],\
                                        const std::vector<struct segment_node*> feasible_segments,\
                                        const std::vector<struct face> obstacle_faces,\
                                        const std::vector<struct pmj_lat> lat_points);
extern "C" struct segment_node* minimize_weighted_activation_time_with_angle_restriction (struct cco_network *the_network,\
                                        struct cost_function_config *config,\
                                        struct local_optimization_config *local_opt_config,\
                                        const double new_pos[],\
                                        const std::vector<struct segment_node*> feasible_segments,\
                                        const std::vector<struct face> obstacle_faces,\
                                        const std::vector<struct pmj_lat> lat_points);
extern "C" struct segment_node* minimize_activation_time_with_angle_restriction_and_using_lat (struct cco_network *the_network,\
                                        struct cost_function_config *config,\
                                        struct local_optimization_config *local_opt_config,\
                                        const double new_pos[],\
                                        const std::vector<struct segment_node*> feasible_segments,\
                                        const std::vector<struct face> obstacle_faces,\
                                        const std::vector<struct pmj_lat> lat_points);
//extern "C" struct segment_node* minimize_weighted_activation_time_with_angle_restriction_and_using_lat (struct cco_network *the_network,\
                                        struct cost_function_config *config,\
                                        struct local_optimization_config *local_opt_config,\
                                        const double new_pos[],\
                                        const std::vector<struct segment_node*> feasible_segments,\
                                        const std::vector<struct face> obstacle_faces,\
                                        const std::vector<struct pmj_lat> lat_points);

// Auxiliary functions

#endif