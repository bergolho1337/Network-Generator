//
// Created by bergolho on 26/01/20.
//

#ifndef NETWORK_H
#define NETWORK_H

#include <cstdio>
#include <cstdlib>
#include <cstdbool>
#include <cstdint>
#include <cstring>
#include <vector>

#include "../utils/utils.h"

struct network_point
{
    uint32_t id;
    double x, y, z;
};

struct network_line
{
    uint32_t src;
    uint32_t dest;
};

void read_initial_network_file (const char filename[], std::vector<struct network_point> &points, std::vector<struct network_line> &lines);
void calculate_unitary_vector(const double x1, const double y1, const double z1, const double x2, const double y2, const double z2,\
                                double u[], double &norm);

#endif