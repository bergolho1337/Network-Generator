#ifndef BRANCH_H
#define BRANCH_H

#include <cstdio>
#include <cstdlib>
#include <cstdbool>
#include <cstdint>
#include <cmath>

#include "../utils/utils.h"

class Branch
{
public:
    uint32_t id;
    double pos[3];
    double dir[3];
    double original_dir[3];
    double length;
    uint32_t counter;
    uint32_t parent;
public:
    Branch (const uint32_t id, const double x, const double y, const double z,\
            const double dx, const double dy, const double dz,\
            uint32_t parent);
    void get_next_branch_position (double new_pos[]);
    void reset ();
    void print ();
};

#endif
