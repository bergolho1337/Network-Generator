#ifndef COST_FUNCTION_H
#define COST_FUNCTION_H

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <string>
#include <vector>

class CCO_Network;
class CostFunctionConfig;
class LocalOptimizationConfig;
class Point;
class Segment;

// Abstract class - Interface
class CostFunction
{
public:
    virtual void init_parameters (CostFunctionConfig *cost_function_config) = 0;
    virtual bool check_restrictions (CCO_Network *the_network, Segment *iconn, Segment *ibiff, Segment *inew) = 0;
    virtual Segment* eval (CCO_Network *the_network,\
                        CostFunctionConfig *cost_function_config,\
                        LocalOptimizationConfig *local_opt_config,\
                        std::vector<Segment*> feasible_segments,\
                        Point *new_term) = 0;
};

#endif