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
#define PRINT_LINE "============================================================================"

static const double ETA = 3.6;          // Blood viscosity
static const double GAMMA = 3.0;        // Bifurcation expoent
static const uint32_t NTOSS = 200;      // Number of tosses for a new terminal
// =================================================================

struct cco_network
{
    bool using_cloud_points;
    char *cloud_points_filename;

    int num_terminals;

    int N_term;
    double Q_perf;
    double p_perf;
    double p_term;
    double r_perf;
    double r_supp;

    double A_perf;

    struct point_list *point_list;
    struct segment_list *segment_list;

    FILE *log_file;
};

void usage (const char pname[]);

struct cco_network* new_cco_network (struct user_options *options);
void free_cco_network (struct cco_network *the_network);

void build_segment (struct cco_network *the_network, const uint32_t index, const double new_pos[]);

void rescale_root (struct segment_node *iroot, const double Q_perf, const double delta_p);
void rescale_tree (struct segment_node *ibiff, const double Q_perf, const double delta_p, const int num_terminals);
void rescale_until_root (struct segment_node *ipar, struct segment_node *ipar_left, struct segment_node *ipar_right,\
                        const double Q_perf, const double delta_p, const int num_terminals);

void recalculate_radius (struct cco_network *the_network);

void calc_middle_point_segment (struct segment_node *s, double pos[]);
void calc_relative_resistance_term (struct segment_node *iterm);
void calc_relative_resistance_subtree (struct segment_node *ibiff, struct segment_node *iconn, struct segment_node *inew);
void calc_radius_term (struct segment_node *iterm, const double Q_term, const double delta_p);
double calc_bifurcation_ratio (const double radius_ratio, bool sign);
double calc_radius_ratio (struct segment_node *iconn, struct segment_node *inew, const double Q_term);
double calc_radius (struct cco_network *the_network, struct segment_node *s);

void check_bifurcation_rule (struct cco_network *the_network);
bool check_collisions (struct cco_network *the_network, const double new_pos[]);

bool has_collision (struct segment_list *s_list, struct segment_node *s, const double new_pos[], FILE *log_file);
bool connection_search (struct cco_network *the_network, const double pos[], const double d_threash);
bool distance_criterion (struct segment_node *s, const double pos[], const double d_threash);

struct segment_node* find_closest_segment (struct cco_network *the_network, const double new_pos[]); 

void make_root (struct cco_network *the_network);
void grow_tree (struct cco_network *the_network);
void generate_terminal (struct cco_network *the_network);
void write_to_vtk (struct cco_network *the_network);

void make_root_using_cloud_points (struct cco_network *the_network, std::vector<struct point> cloud_points);
uint32_t sort_point_from_cloud (double pos[], std::vector<struct point> cloud_points);
void read_cloud_points (const char filename[], std::vector<struct point> &cloud_points);
void generate_terminal_using_cloud_points(struct cco_network *the_network, std::vector<struct point> cloud_points);


// Test functions
void test1 (struct cco_network *the_network);
void test2 (struct cco_network *the_network); 
void test3 (struct cco_network *the_network);
void test_cco (struct cco_network *the_network);
void test_cco_using_cloud (struct cco_network *the_network);

#endif