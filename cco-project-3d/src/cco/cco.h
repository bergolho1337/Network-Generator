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

#include "../point-list/point-list.h"
#include "../segment-list/segment-list.h"
#include "../options/user_options.h"
#include "../utils/utils.h"

#include "cco_helper.h"

// CONSTANTS AND MACROS 
// =================================================================
static const double ETA = 3.6e-03;                  // Blood viscosity 
static const double GAMMA = 3.0;                    // Bifurcation expoent
static const uint32_t NTOSS = 10;                   // Number of tosses for a new terminal
static const double FACTOR = 0.95;                  // Reduction factor for the distance criterion

#define PRINT_LINE "========================================================================================================"
// =================================================================

struct cco_network
{
    int num_terminals;                      

    int N_term;                         // Total number of terminals 
    double Q_perf;                      // Root perfusion
    double p_perf;                      // Perfusion pressure over the root
    double p_term;                      // Perfusion pressure over the terminals
    double r_perf;                      // Perfusion radius
    double r_supp;                                            

    double A_perf;
    double V_perf;                      // Perfusion volume (sphere)

    double Q_term;                      // Terminals perfusion
    double delta_p;                     // Pressure drop

    double root_pos[3];

    struct point_list *point_list;
    struct segment_list *segment_list;

    bool using_cloud_points;
    char *cloud_points_filename;

    bool using_local_optimization;
    char *local_optimization_function_name;

    char *cost_function_name;

    FILE *log_file;
};

void usage (const char pname[]);

struct cco_network* new_cco_network (struct user_options *options);
void free_cco_network (struct cco_network *the_network);

void set_parameters (struct cco_network *the_network, struct user_options *options);
void set_cost_function_name (struct cco_network *the_network, struct user_options *options);
void set_cloud_points_name (struct cco_network *the_network, struct user_options *options);
void set_local_optimization_function_name (struct cco_network *the_network, struct user_options *options);

struct segment_node* build_segment (struct cco_network *the_network, struct local_optimization_config *local_opt_config,\
                                const uint32_t index, const double new_pos[]);

void rescale_root (struct segment_node *iroot, const double Q_perf, const double delta_p);
void rescale_tree (struct segment_node *ibiff, struct segment_node *iconn, struct segment_node *inew,\
                 const double Q_perf, const double delta_p, const int num_terminals);
void rescale_until_root (struct segment_node *ipar, struct segment_node *ipar_left, struct segment_node *ipar_right,\
                        const double Q_perf, const double delta_p, const int num_terminals);

void recalculate_radius (struct cco_network *the_network);

void restore_state_tree (struct cco_network *the_network,\
                        struct segment_node *iconn);

void check_bifurcation_rule (struct cco_network *the_network);
bool check_collisions_and_fill_feasible_segments (struct cco_network *the_network, const double new_pos[],\
                    std::vector<struct segment_node*> &feasible_segments);

bool has_collision (struct segment_list *s_list, struct segment_node *s, const double new_pos[], FILE *log_file);
bool has_collision (struct segment_list *s_list, struct segment_node *iconn, struct segment_node *ibiff, struct segment_node *inew, FILE *log_file);

bool connection_search (struct cco_network *the_network, const double pos[], const double d_threash);
bool distance_criterion (struct segment_node *s, const double pos[], const double d_threash);

void grow_tree (struct cco_network *the_network, struct user_options *options);
void grow_tree_using_cloud_points (struct cco_network *the_network, struct user_options *options, std::vector<struct point> cloud_points);

void make_root_using_cloud_points (struct cco_network *the_network, std::vector<struct point> cloud_points);
                        
void generate_terminal_using_cloud_points(struct cco_network *the_network,\
                                          struct cost_function_config *config,\
                                          struct local_optimization_config *local_opt_config,\
                                          std::vector<struct point> cloud_points);

uint32_t sort_point_from_cloud (double pos[], std::vector<struct point> cloud_points);
void sort_point_from_cloud_v2 (double pos[], std::vector<struct point> cloud_points);

void read_cloud_points (const char filename[], std::vector<struct point> &cloud_points);
void build_cloud_points (std::vector<struct point> &cloud_points, const double radius);

#endif