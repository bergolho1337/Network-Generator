#ifndef ACTIVATION_TIME_FUNCTION_H
#define ACTIVATION_TIME_FUNCTION_H

#include "cost_function.h"

class ActivationTimeFunction : public CostFunction
{
public:
    void init_parameters (CostFunctionConfig *cost_function_config);
    bool check_restrictions (CCO_Network *the_network, Segment *iconn, Segment *ibiff, Segment *inew);
    Segment* eval (CCO_Network *the_network,\
                        CostFunctionConfig *cost_function_config,\
                        LocalOptimizationConfig *local_opt_config,\
                        std::vector<Segment*> feasible_segments,\
                        Point *new_term);
};

#endif