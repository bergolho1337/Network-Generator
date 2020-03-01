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

#include <vtkRegularPolygonSource.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkSphereSource.h>

#include "../cco/cco.h"

#include "dSFMT.h"

// =================================================================
// CONSTANTS AND MACROS

#define UM_TO_CM 0.0001
#define CM_TO_UM 10000.0
#define MM_TO_UM 1000.0
#define MS_TO_S 0.001
#define S_TO_MS 1000.0
#define CM_S_TO_M_S 0.01
#define MS_TO_US 1000.0
#define CM_TO_MM 10.0

#define PRINT_DOTS ".........................................................................."
#define PRINT_STARS "******************************************************************************"

static const double EPSILON = 1.0e-02;                   // Tolerance for comparing real numbers
static const uint32_t RAND_ARRAY_SIZE = 8000000;         // Size of randomic vector
static const uint32_t RAND_SEED = 1;                     // Random seed
static const uint32_t TOTAL_CLOUD_POINTS_SIZE = 2500000; // Total number of points in the generated cloud of points
// =================================================================

struct random_generator
{
    uint32_t counter;
    double *array;
};

struct random_generator* new_random_generator ();
void free_random_generator (struct random_generator *the_generator);
void generate_random_array (struct random_generator *the_generator);
double get_value (struct random_generator *the_generator);

double generate_random_number_default ();

void generate_cloud_points (struct random_generator *the_generator, std::vector<struct point> &cloud_points, const double radius);


bool collision_detection (const double x1, const double y1, const double z1,\
                          const double x2, const double y2, const double z2,\
                          const double x3, const double y3, const double z3,\
                          const double x4, const double y4, const double z4);

void build_unitary_vector (double u[], const double x1, const double y1, const double z1,\
                                       const double x2, const double y2, const double z2);

double calc_angle_between_vectors (const double u[], const double v[]);

double euclidean_norm (const double x1, const double y1, const double z1,\
                    const double x2, const double y2, const double z2);
double calc_dot_product (const double u[], const double v[]);

double calc_dthreashold (const double radius, const int num_terminals);
double calc_dproj (struct segment_node *s, const double pos[]);
double calc_dortho (struct segment_node *s, const double pos[]);
double calc_dend (struct segment_node *s, const double pos[]);

double calc_perfusion_radius (const double V);
double calc_flux_terminals (const double Q_perf, const int N_term);

void calc_plane_coefficients (const double v1[], const double v2[], const double v3[], double N[], double &D);
void calc_subtract_vector(const double x_new[], const double x_prox[], double seg_rq[]);

void draw_perfusion_volume (const double radius);

bool check_size (const double p[]);

bool check_segment_plane_intersection (const double x_prox[], const double x_new[], struct face the_face);

void print_terminal_activation_time (struct cco_network *the_network,\
                            const double c, const double cm, const double rc, const double rm);

void write_to_vtk (struct cco_network *the_network);
void write_to_vtk_iteration (struct cco_network *the_network);

#endif
