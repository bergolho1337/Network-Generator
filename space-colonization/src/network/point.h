#ifndef POINT_H
#define POINT_H

#include <cstdio>
#include <cstdlib>
#include <cstdbool>
#include <cstdint>
#include <cmath>

#include "../utils/utils.h"

class Point
{
public:
    uint32_t index;
    double pos[3];
public:
    Point (const uint32_t index, const double x, const double y, const double z);
    void print ();
};

#endif
