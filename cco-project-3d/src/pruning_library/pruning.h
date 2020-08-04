#ifndef PRUNING_LIB_H
#define PRUNING_LIB_H

#include <cstdio>
#include <cstdlib>
#include <cmath>

#include "../options/pruning_config.h"
#include "../utils/utils.h"

// Pruning functions headers
extern "C" double hyperbolic_tangent (struct pruning_config *config, const double level, const double length);
extern "C" double exponential (struct pruning_config *config, const double level, const double length);
extern "C" double length (struct pruning_config *config, const double level, const double length);
extern "C" double hyperbolic_tangent_with_length (struct pruning_config *config, const double level, const double length);

// Auxiliary functions


#endif