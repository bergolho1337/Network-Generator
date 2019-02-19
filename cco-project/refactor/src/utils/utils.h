//
// Created by bergolho on 19/02/19.
//

#ifndef UTILS_H
#define UTILS_H

#include <cstdio>
#include <cstdlib>
#include <cstdbool>
#include <cstdint>
#include <cmath>

void generate_point_inside_perfusion_area (double pos[], const double radius);

double generate_random_number ();

double euclidean_norm (const double x1, const double y1, const double z1,\
                    const double x2, const double y2, const double z2);

#endif