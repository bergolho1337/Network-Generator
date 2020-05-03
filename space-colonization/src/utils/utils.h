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
inline double generate_random_number () { return (double)rand() / (double)RAND_MAX; }
inline void calculate_difference_vector (const double a[], const double b[], double c[]) { for (uint32_t i = 0; i < 3; i++) c[i] = a[i] - b[i]; }
inline double calculate_distance (const double x1, const double y1, const double z1,\
                           const double x2, const double y2, const double z2) { return sqrt( powf(x2-x1,2) + powf(y2-y1,2) + powf(z2-z1,2) ); }
inline void calculate_unitary_vector (const double a[], const double b[], double u[]) 
{ 
    double norm = calculate_distance(a[0],a[1],a[2],b[0],b[1],b[2]);
    for (uint32_t i = 0; i < 3; i++)
        u[i] = (b[i] - a[i]) / norm; 
}
inline void normalize_vector (double a[])  
{ 
    double norm = sqrt( pow(a[0],2) + pow(a[1],2) + pow(a[2],2));
    for (uint32_t i = 0; i < 3; i++)
        a[i] /= norm;
}

void usage (const char pname[]);

#endif
