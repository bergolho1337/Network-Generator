#ifndef SQUARE_LIB_H
#define SQUARE_LIB_H

#include <cstdio>
#include <cstdlib>
#include <cmath>

#include "../generator/generator.h"

// Cost functions implementations
extern "C" void default_cloud_square (struct cloud_generator_data *generator,\
                                      struct user_data *config);
extern "C" void spaced_cloud_square (struct cloud_generator_data *generator,\
                                     struct user_data *config);


// Auxiliary functions
void draw_square_area (const double side_length);
void generate_point_inside_square (double pos[], const double side_length);

#endif