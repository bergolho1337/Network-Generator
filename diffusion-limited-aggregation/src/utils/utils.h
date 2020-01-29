
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

void usage (const char pname[]);
double calculate_distance (const double p1[], const double p2[]);
double generate_random_number ();
bool get_parameter_value_from_map (std::map<std::string,double> *params,\
                                const std::string key, double *value);
void print_progress_bar (const uint32_t cur_iter, const uint32_t max_iter);

#endif