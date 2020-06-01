#ifndef LEAF_H
#define LEAF_H

#include <cstdio>
#include <cstdlib>
#include <cstdbool>
#include <cstdint>
#include <cmath>

#include "../utils/utils.h"

class Leaf
{
public:
    bool is_reached;
    double pos[3];
public:
    Leaf (const double x, const double y, const double z);
    void print ();
};

#endif
