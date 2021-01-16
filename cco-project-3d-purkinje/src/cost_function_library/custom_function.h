#ifndef CUSTOM_FUNCTION_H
#define CUSTOM_FUNCTION_H

#include "cost_function.h"

class MinimizeCustomFunction : public CostFunction
{
public:
    double beta;
    double alpha;
    double min_degrees_limit;
    double max_degrees_limit;
    double min_segment_length;
    double max_segment_length;
public:
    void init_parameters (CostFunctionConfig *cost_function_config);
    bool check_restrictions (CCO_Network *the_network, Segment *iconn, Segment *ibiff, Segment *inew);
    bool check_restrictions_2 (CCO_Network *the_network, Segment *iconn, Segment *ibiff, Segment *inew);
    Segment* eval (CCO_Network *the_network,\
                        CostFunctionConfig *cost_function_config,\
                        LocalOptimizationConfig *local_opt_config,\
                        std::vector<Segment*> feasible_segments,\
                        Point *new_term, const bool is_pmj);
    double calc_custom_function (CCO_Network *the_network);
    double calc_segment_custom_function (Segment *s);
    bool check_angle_restriction (Segment *iconn, Segment *inew);
    bool check_minimum_segment_size (Segment *iconn, Segment *ibiff, Segment *inew);
    bool check_maximum_segment_size (Segment *iconn, Segment *ibiff, Segment *inew);
    bool check_collision (CCO_Network *the_network, Segment *iconn, Segment *ibiff, Segment *inew);
};

class MaximizeCustomFunction : public CostFunction
{
public:
    double beta;
    double alpha;
    double min_degrees_limit;
    double max_degrees_limit;
    double min_segment_length;
    double max_segment_length;
public:
    void init_parameters (CostFunctionConfig *cost_function_config);
    bool check_restrictions (CCO_Network *the_network, Segment *iconn, Segment *ibiff, Segment *inew);
    bool check_restrictions_2 (CCO_Network *the_network, Segment *iconn, Segment *ibiff, Segment *inew);
    Segment* eval (CCO_Network *the_network,\
                        CostFunctionConfig *cost_function_config,\
                        LocalOptimizationConfig *local_opt_config,\
                        std::vector<Segment*> feasible_segments,\
                        Point *new_term, const bool is_pmj);
    double calc_custom_function (CCO_Network *the_network);
    double calc_segment_custom_function (Segment *s);
    bool check_angle_restriction (Segment *iconn, Segment *inew);
    bool check_minimum_segment_size (Segment *iconn, Segment *ibiff, Segment *inew);
    bool check_maximum_segment_size (Segment *iconn, Segment *ibiff, Segment *inew);
    bool check_collision (CCO_Network *the_network, Segment *iconn, Segment *ibiff, Segment *inew);
};

#endif