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
class Segment;
class Point;

void create_directory (const char *path);
void calc_mean_std (std::vector<double> arr, double &mean, double &std);
void write_vector_to_file (std::vector<double> arr, std::string filename);
void write_geometric_info_to_file (std::string filename,\
                    std::vector<double> segments,const double mean_segment_length, const double std_segment_length,\
                    std::vector<double> angles,const double mean_biff_angle, const double std_biff_angle);
void write_electric_info_to_file (std::string filename,\
                    const double max_lat_error, const double min_ref_lat, const double max_ref_lat,\
                    const double min_aprox_lat, const double max_aprox_lat,\
                    const double rmse, const double rrmse,\
                    const double epsilon_2ms, const double epsilon_5ms);

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

#endif
