#ifndef UTILS_H_
#define UTILS_H_

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <cstdint>

#define PRINT_LINE "---------------------------------------------------------------------------------"

void usage (const char pname[]);
void print_stars (const int number);
double calc_angle_between_vectors (const double u[], const double v[]);
double calc_norm (double x1, double y1, double z1, double x2, double y2, double z2);
bool check_file_extension (const char filename[], const char extension_name[]);


#endif