#include "pruning.h"

double pruning_function (const double level, const double A, const double B, const double C)
{
    return 50.0*tanh(-0.25*level+3.0) + 50.0;

    //else
    //    return A*expf(B*level) + C;
}