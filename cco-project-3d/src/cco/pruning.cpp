#include "pruning.h"

double pruning_function (const double level, const double A, const double B, const double C)
{
    if (level < 5)
        return 100.0;
    else
        return A*expf(B*level) + C;
}