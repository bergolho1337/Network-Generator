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
#include <sys/stat.h> 
#include <sys/types.h> 

#include <string>
#include <sstream>

#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkUnstructuredGrid.h>
#include <vtkCell.h>
#include <vtkLine.h>
#include <vtkHexahedron.h>
#include <vtkVertex.h>
#include <vtkFloatArray.h>
#include <vtkCellData.h>
#include <vtkPointData.h>
#include <vtkGenericDataObjectReader.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkDijkstraGraphGeodesicPath.h>
#include <vtkSTLReader.h>

#include "../cco/cco.h"

#include "dSFMT.h"

// =================================================================
// CONSTANTS AND MACROS
#define UM_TO_CM 0.0001
#define CM_TO_UM 1.0E+05
#define MM_TO_UM 1.0E+03
#define M_TO_UM 1.0E+06
#define M_TO_MM 1.0E+03
#define MS_TO_S 0.001
#define S_TO_MS 1000.0
#define CM_S_TO_M_S 0.01
#define MS_TO_US 1000.0
#define CM_TO_MM 10.0
#define CM_TO_M 0.01

static const double EPSILON = 1.0e-02;                   // Tolerance for comparing real numbers
static const uint32_t RAND_ARRAY_SIZE = 8000000;         // Size of randomic vector
static const uint32_t RAND_SEED = 1;                     // Random seed
static const uint32_t TOTAL_CLOUD_POINTS_SIZE = 2500000; // Total number of points in the generated cloud of points
// =================================================================

class CCO_Network;

void create_directory (const char *path);
void calc_mean_std (std::vector<double> arr, double &mean, double &std);
void write_vector_to_file (std::vector<double> arr, std::string filename);
void write_info_to_file (std::string filename,\
                    std::vector<double> segments,const double mean_segment_length, const double std_segment_length,\
                    std::vector<double> angles,const double mean_biff_angle, const double std_biff_angle);

double euclidean_norm (const double x1, const double y1, const double z1,\
                    const double x2, const double y2, const double z2);
double calc_angle_between_vectors (const double u[], const double v[]);

double calc_dthreashold (const double radius, const uint32_t num_terminals);
double calc_dproj (Segment *s, const double pos[]);
double calc_dortho (Segment *s, const double pos[]);
double calc_dend (Segment *s, const double pos[]);
double calc_dot_product (const double u[], const double v[]);

bool read_points_from_vtk (const char filename[], std::vector<Point*> &points);

bool check_size (const double p[]);
bool collision_detection (const double x1, const double y1, const double z1,\
                          const double x2, const double y2, const double z2,\
                          const double x3, const double y3, const double z3,\
                          const double x4, const double y4, const double z4);
/*
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




void build_unitary_vector (double u[], const double x1, const double y1, const double z1,\
                                       const double x2, const double y2, const double z2);






double calc_perfusion_radius (const double V);
double calc_flux_terminals (const double Q_perf, const int N_term);

void calc_plane_coefficients (const double v1[], const double v2[], const double v3[], double N[], double &D);
void calc_subtract_vector(const double x_new[], const double x_prox[], double seg_rq[]);

void draw_perfusion_volume (const double radius);

bool check_segment_plane_intersection (const double x_prox[], const double x_new[], struct face the_face);

void print_terminal_activation_time (struct cco_network *the_network,\
                            const double G, const double Cf, const double tau_f);
                            
*/

#endif
