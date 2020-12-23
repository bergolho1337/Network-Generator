#include "activation_time_function.h"

#include "../cco/cco.h"

void ActivationTimeFunction::init_parameters (CostFunctionConfig *cost_function_config)
{
    // Parameters of the function goes here ...
}

bool ActivationTimeFunction::check_restrictions (CCO_Network *the_network, Segment *iconn, Segment *ibiff, Segment *inew)
{
    return true;    
}

Segment* ActivationTimeFunction::eval (CCO_Network *the_network,\
                        CostFunctionConfig *cost_function_config,\
                        LocalOptimizationConfig *local_opt_config,\
                        std::vector<Segment*> feasible_segments,\
                        Point *new_term)
{
    printf("Activation time function\n");
    return NULL;
}