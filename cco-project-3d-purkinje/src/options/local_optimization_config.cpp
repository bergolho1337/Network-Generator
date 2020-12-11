#include "local_optimization_config.h"

#include "../segment/segment.h"
#include "../point/point.h"

LocalOptimizationConfig::LocalOptimizationConfig ()
{
    this->first_call = true;
}

LocalOptimizationConfig::~LocalOptimizationConfig ()
{

}

void LocalOptimizationConfig::save_original_bifurcation_position (Segment *s, double ori_pos[])
{
    ori_pos[0] = s->src->x;
    ori_pos[1] = s->src->y;
    ori_pos[2] = s->src->z;
}

void LocalOptimizationConfig::update_bifurcation_position (Point *p)
{
    this->best_pos[0] = p->x;
    this->best_pos[1] = p->y;
    this->best_pos[2] = p->z;
}

void LocalOptimizationConfig::initialize_best_position_as_middle_point (double best_pos[], const double ori_pos[])
{
    best_pos[0] = ori_pos[0];
    best_pos[1] = ori_pos[1];
    best_pos[2] = ori_pos[2];
}

void LocalOptimizationConfig::print ()
{
    printf("Local optimization function name = \"%s\"\n",this->function_name.c_str());
}