#ifndef PRUNING_LIB_H
#define PRUNING_LIB_H

#include <cstdio>
#include <cstdlib>
#include <cmath>

#include "../options/pruning_config.h"
#include "../utils/utils.h"

// Pruning functions headers
extern "C" double hyperbolic_tangent (struct pruning_config *config, const double level);
extern "C" double exponential (struct pruning_config *config, const double level);

// Auxiliary functions


#endif