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

// =================================================================
// CONSTANTS AND MACROS
#define PRINT_LINE "--------------------------------------------------------------------------------"
#define PRINT_DOTS "................................................................................"
#define PRINT_STARS "*******************************************************************************"

inline double convert_degrees_to_radians (const double degrees) { return degrees * M_PI / 180.0; }
inline double calculate_distance (const double x1, const double y1, const double z1,\
                           const double x2, const double y2, const double z2) { return sqrt( powf(x2-x1,2) + powf(y2-y1,2) + powf(z2-z1,2) ); }

void usage (const char pname[]);

void calculate_unitary_vector (const double a[], const double b[], double u[]);
void calculate_new_branch_position (double new_pos[], const double cur_pos[], const double u[], const double length, const double angle_degree);

#endif
