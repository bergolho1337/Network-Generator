
#ifndef UTILS_H
#define UTILS_H

#include <cstdio>
#include <cstdlib>
#include <cstdbool>
#include <cstdint>
#include <cmath>

#include <vector>
#include <map>
#include <string>

#define PRINT_LINE "==========================================================================================="
#define PRINT_DOTS "..........................................................................................."
#define EPSILON 1.0E-08

#define OI printf("oi\n")

void usage (const char pname[]);
double calculate_distance (const double p1[], const double p2[]);
double calculate_euclidean_norm (const double x1, const double y1, const double z1,\
                                const double x2, const double y2, const double z2);
double calculate_angle_between_vectors (const double u[], const double v[]);

double generate_random_number ();
bool get_parameter_value_from_map (std::map<std::string,double> *params,\
                                const std::string key, double *value);
void print_progress_bar (const uint32_t cur_iter, const uint32_t max_iter);

#endif