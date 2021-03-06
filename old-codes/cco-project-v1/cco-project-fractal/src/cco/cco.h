//
// Created by bergolho on 12/02/19.
//

#ifndef CCO_H
#define CCO_H

#include <cstdio>
#include <cstdlib>
#include <cstdbool>
#include <cstdint>
#include <cmath>

#include <vector>
#include <algorithm>

#include <vtkRegularPolygonSource.h>
#include <vtkXMLPolyDataWriter.h>

#include "../point-list/point-list.h"
#include "../segment-list/segment-list.h"
#include "../options/user_options.h"
#include "../utils/utils.h"

// CONSTANTS AND MACROS 
// =================================================================
static const double ETA = 3.6;          // Blood viscosity (cP)
static const double GAMMA = 3.0;        // Bifurcation expoent
static const uint32_t NTOSS = 200;      // Number of tosses for a new terminal
static const double MEAN = 1.0;         // Mean value for the normal distribution
static const double STDDEV = 0.0;       // Standart deviation for the normal distribution

#define PRINT_LINE "========================================================================================================"
// =================================================================

struct cco_network
{
    int num_terminals;

    int N_term;
    double start_radius;

    double Q_perf;
    double p_perf;
    double p_term;
    double r_perf;
    double r_supp;

    double A_perf;

    struct point_list *point_list;
    struct segment_list *segment_list;

    bool using_cloud_points;
    char *cloud_points_filename;

    char *cost_function_name;

    FILE *log_file;
};

void usage (const char pname[]);

struct cco_network* new_cco_network (struct user_options *options);
void free_cco_network (struct cco_network *the_network);

struct segment_node* build_segment (struct cco_network *the_network, const uint32_t index, const double new_pos[]);

void rescale_root (struct segment_node *iroot, const double Q_perf, const double delta_p, const double start_radius);
void rescale_tree (struct segment_node *ibiff, struct segment_node *iconn, struct segment_node *inew,\
                 const double Q_perf, const double delta_p, const double start_radius, const int num_terminals);
void rescale_until_root (struct segment_node *ipar, struct segment_node *ipar_left, struct segment_node *ipar_right,\
                        const double Q_perf, const double delta_p, const double start_radius, const int num_terminals);

void recalculate_radius (struct cco_network *the_network);

void restore_state_tree (struct cco_network *the_network,\
                        struct segment_node *iconn);

void calc_middle_point_segment (struct segment_node *s, double pos[]);
void calc_relative_resistance_term (struct segment_node *iterm);
void calc_relative_resistance_subtree (struct segment_node *ibiff, struct segment_node *iconn, struct segment_node *inew);
void calc_pressure_drop_term (struct segment_node *iterm, const double Q_term);
void calc_pressure_drop_subtree (struct segment_node *iconn, const double Q_term);
void calc_radius_term (struct segment_node *iterm, const double Q_term, const double delta_p);
double calc_bifurcation_ratio (const double radius_ratio, bool sign);
double calc_radius_ratio (struct segment_node *iconn, struct segment_node *inew, const double Q_term);
double calc_radius (struct cco_network *the_network, struct segment_node *s);

double calc_tree_volume (struct cco_network *the_network);
double calc_segment_volume (struct segment_node *s);
double calc_assymetric_ratio (struct segment_node *right, struct segment_node *left);

double calc_terminal_activation_time (struct segment_node *s,\
                        const double c, const double cm, const double rc, const double rm);
double calc_segment_activation_time (struct segment_node *s,\
                        const double c, const double cm, const double rc, const double rm);
double calc_propagation_velocity (const double r,\
                        const double c, const double cm, const double rc, const double rm);
double calc_lambda_m (const double r, const double rc, const double rm);
double calc_tau_m (const double cm, const double rm);

double calc_segment_activation_time_using_level (const double at, struct segment_node *iconn);
double calc_segment_level (struct segment_node *iconn);

double calc_segment_custom_function (struct segment_node *s, const double beta, const double alpha);
double calc_custom_function (struct cco_network *the_network, const double beta, const double alpha);

void check_bifurcation_rule (struct cco_network *the_network);
bool check_collisions (struct cco_network *the_network, const double new_pos[],\
                    std::vector<struct segment_node*> &feasible_segments);

bool has_collision (struct segment_list *s_list, struct segment_node *s, const double new_pos[], FILE *log_file);
bool connection_search (struct cco_network *the_network, const double pos[], const double d_threash);
bool distance_criterion (struct segment_node *s, const double pos[], const double d_threash);

bool has_deviation (struct segment_list *s_list, struct segment_node *inew,\
                    const double new_at, const double limit,\
                    const double c, const double cm, const double rc, const double rm);
bool is_terminal (struct segment_node *s);

void grow_tree (struct cco_network *the_network, struct user_options *options);
void grow_tree_default (struct cco_network *the_network, struct user_options *options);
void grow_tree_using_cloud_points (struct cco_network *the_network, struct user_options *options);

void make_root (struct cco_network *the_network);
void make_root_using_cloud_points (struct cco_network *the_network, std::vector<struct point> cloud_points);

void generate_terminal (struct cco_network *the_network, struct cost_function_config *config);

uint32_t sort_point_from_cloud (double pos[], std::vector<struct point> cloud_points);
void read_cloud_points (const char filename[], std::vector<struct point> &cloud_points);
void generate_terminal_using_cloud_points(struct cco_network *the_network, struct cost_function_config *config,\
                                         std::vector<struct point> cloud_points);

void write_to_vtk (struct cco_network *the_network);

#endif