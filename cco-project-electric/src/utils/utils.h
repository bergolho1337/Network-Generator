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

#include "../cco/cco.h"

// =================================================================
// CONSTANTS AND MACROS 

#define UM_TO_CM 0.0001
#define CM_TO_UM 10000.0 
#define MS_TO_S 0.001
#define S_TO_MS 1000.0
#define CM_S_TO_M_S 0.01
#define MS_TO_US 1000.0

#define PRINT_DOTS ".........................................................................."
#define PRINT_STARS "******************************************************************************"

// =================================================================

void generate_point_inside_perfusion_area (double pos[], const double radius);

double generate_random_number ();

bool collision_detection (const double x1, const double y1, const double z1,\
                          const double x2, const double y2, const double z2,\
                          const double x3, const double y3, const double z3,\
                          const double x4, const double y4, const double z4);

void build_unitary_vector (double u[], const double x1, const double y1, const double z1,\
                                       const double x2, const double y2, const double z2);

double calc_angle_between_vectors (const double u[], const double v[]);

double euclidean_norm (const double x1, const double y1, const double z1,\
                    const double x2, const double y2, const double z2);

double calc_dthreashold (const double radius, const int num_terminals);
double calc_dproj (struct segment_node *s, const double pos[]);
double calc_dortho (struct segment_node *s, const double pos[]);
double calc_dend (struct segment_node *s, const double pos[]);

void draw_perfusion_area (struct cco_network *the_network);

void print_terminal_activation_time (struct cco_network *the_network,\
                            const double c, const double cm, const double rc, const double rm);

void write_to_vtk (struct cco_network *the_network);

#endif