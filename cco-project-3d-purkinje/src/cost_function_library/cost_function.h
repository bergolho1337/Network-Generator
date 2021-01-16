#ifndef COST_FUNCTION_H
#define COST_FUNCTION_H

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
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
    virtual bool check_restrictions_2 (CCO_Network *the_network, Segment *iconn, Segment *ibiff, Segment *inew) = 0;
    virtual Segment* eval (CCO_Network *the_network,\
                        CostFunctionConfig *cost_function_config,\
                        LocalOptimizationConfig *local_opt_config,\
                        std::vector<Segment*> feasible_segments,\
                        Point *new_term, const bool is_pmj) = 0;
};

class Evaluation
{
public:
    double eval;
    double error;
    double best_pos[3];
    Segment *best;
public:
    Evaluation (const double eval, const double error, const double best_pos[], Segment *best)
    {
        this->error = error;
        this->eval = eval;
        memcpy(this->best_pos,best_pos,sizeof(double)*3);
        this->best = best;
    }
};

#endif