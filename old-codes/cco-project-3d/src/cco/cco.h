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
#include "../face-list/face-list.h"
#include "../options/user_options.h"
#include "../utils/utils.h"
#include "../utils/stop_watch.h"

#include "cco_helper.h"

// CONSTANTS AND MACROS 
// =================================================================
static const double ETA = 3.6e-03;                  // Blood viscosity 
static const uint32_t NTOSS = 10;                   // Number of tosses for a new terminal
static const double FACTOR = 0.95;                  // Reduction factor for the distance criterion
static const double EPISON = 1.0e-16;                // Minimum size for a segment      CHANGE
static const uint32_t PRUNING_PASSES = 1;           // Number of times the pruning procedure will be called
#define PRINT_LINE "========================================================================================================"
// =================================================================

struct cco_network
{
    uint32_t seed;                      // Random seed
    uint32_t max_rand_offset;           // Maximum offset for the random generator

    int num_terminals;                      

    double gamma;                       // Bifurcation expoent

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
    bool using_obstacle;
    char *obstacle_filename;
    bool using_pmj_location;
    char *pmj_location_filename;
    bool using_lat;
    char *lat_filename;

    bool using_local_optimization;
    char *local_optimization_function_name;

    bool using_only_murray_law;

    char *cost_function_name;

    char *output_dir;

    bool using_pruning;
    char *pruning_function_name;

    bool using_initial_network;
    char *initial_network_filename;

    FILE *log_file;
};

void usage (const char pname[]);

struct cco_network* new_cco_network (struct user_options *options);
void free_cco_network (struct cco_network *the_network);

void set_parameters (struct cco_network *the_network, struct user_options *options);
void set_save_network (struct cco_network *the_network, struct user_options *options);
void set_cost_function_name (struct cco_network *the_network, struct user_options *options);
void set_cloud_points_name (struct cco_network *the_network, struct user_options *options);
void set_obstacle_name (struct cco_network *the_network, struct user_options *options);
void set_pmj_location_name (struct cco_network *the_network, struct user_options *options);
void set_lat_name (struct cco_network *the_network, struct user_options *options);
void set_local_optimization_function_name (struct cco_network *the_network, struct user_options *options);
void set_pruning_function (struct cco_network *the_network, struct user_options *options);

struct segment_node* build_segment (struct cco_network *the_network, struct local_optimization_config *local_opt_config,\
                                const uint32_t index, const double new_pos[]);

void rescale_root (struct segment_node *iroot, const double Q_perf, const double delta_p, const bool using_only_murray);
void rescale_tree (struct segment_node *ibiff, struct segment_node *iconn, struct segment_node *inew,\
                 const double Q_perf, const double delta_p, const double gamma, const int num_terminals, const bool use_only_murray);
void rescale_until_root (struct segment_node *ipar, struct segment_node *ipar_left, struct segment_node *ipar_right,\
                        const double Q_perf, const double delta_p, const double gamma, const int num_terminals, const bool use_only_murray);

void recalculate_radius (struct cco_network *the_network);

void restore_state_tree (struct cco_network *the_network,\
                        struct segment_node *iconn);

void check_bifurcation_rule (struct cco_network *the_network);
bool check_collisions_and_fill_feasible_segments (struct cco_network *the_network, const double new_pos[],\
                    std::vector<struct segment_node*> &feasible_segments);

bool has_collision (struct segment_list *s_list, struct segment_node *s, const double new_pos[], FILE *log_file);
bool has_collision (struct segment_list *s_list, struct segment_node *iconn, struct segment_node *ibiff, struct segment_node *inew, FILE *log_file);
bool has_intersect_obstacle (struct segment_node *inew, std::vector<struct face> obstacle_faces);
bool has_intersect_obstacle (const double x_prox[], const double x_new[], std::vector<struct face> obstacle_faces);
bool has_valid_segment_sizes (const double iconn_size, const double ibiff_size, const double inew_size);
bool has_valid_segment_sizes_2 (const double iconn_size, const double ibiff_size, const double inew_size);
bool has_level (struct segment_node *inew, const uint32_t num_terminals);

bool connection_search (struct cco_network *the_network, const double pos[], const double d_threash);
bool distance_criterion (struct segment_node *s, const double pos[], const double d_threash);

void grow_tree (struct cco_network *the_network, struct user_options *options);
void grow_tree_using_cloud_points (struct cco_network *the_network,\
                                   struct user_options *options,\
                                   std::vector<struct point> cloud_points,\
                                   std::vector<struct face> obstacle_faces,\
                                   std::vector<struct point> pmj_points,\
                                   std::vector<struct pmj_lat> lat_points);
void generate_terminal_using_pmj_points(struct cco_network *the_network,\
                                          struct cost_function_config *config,\
                                          struct local_optimization_config *local_opt_config,\
                                          std::vector<struct point> pmj_points,\
                                          std::vector<struct face> obstacle_faces,\
                                          std::vector<struct pmj_lat> lat_points);

void make_root_using_cloud_points (struct cco_network *the_network, std::vector<struct point> cloud_points, std::vector<struct face> obstacle_faces);
void make_root_using_initial_network (struct cco_network *the_network);
                        
void generate_terminal_using_cloud_points(struct cco_network *the_network,\
                                          struct cost_function_config *config,\
                                          struct local_optimization_config *local_opt_config,\
                                          std::vector<struct point> cloud_points,\
                                          std::vector<struct face> obstacle_faces,\
                                          std::vector<struct pmj_lat> lat_points);

void prune_tree (struct cco_network *the_network, struct pruning_config *config);
void prune_tree_2 (struct cco_network *the_network, struct pruning_config *config);
void prune_tree_segment (struct cco_network *the_network, struct segment_node *inew);

void sort_point_from_cloud_v1 (double pos[], std::vector<struct point> cloud_points);
void sort_point_from_cloud_v2 (double pos[], std::vector<struct point> cloud_points);
void sort_point_from_cloud_v3 (double pos[], std::vector<struct point> cloud_points);
void sort_point_from_cloud_v4 (double pos[], std::vector<struct point> cloud_points, uint32_t max_rand_offset);

void read_cloud_points (const char filename[], std::vector<struct point> &cloud_points);
void build_cloud_points (std::vector<struct point> &cloud_points, const double radius);

void read_obstacle_faces (const char filename[], std::vector<struct face> &obstacle_faces);
void read_face (FILE *file, std::vector<struct face> &faces);

void read_pmj_location_points (const char filename[], std::vector<struct point> &pmj_points);

void read_lat_points (const char filename[], std::vector<struct pmj_lat> &lat_points);

void print_network_info (struct cco_network *the_network);

#endif