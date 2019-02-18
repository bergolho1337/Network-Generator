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

#include "../point-list/point-list.h"
#include "../segment-list/segment-list.h"
#include "../options/user_options.h"

// CONSTANTS AND MACROS 
// =================================================================
static const double ETA = 3.6;          // Blood viscosity
static const double GAMMA = 3.0;        // Bifurcation expoent
// =================================================================

struct cco_network
{
    int num_terminals;

    int N_term;
    double Q_perf;
    double p_perf;
    double p_term;
    double r_perf;
    double r_supp;

    struct point_list *point_list;
    struct segment_list *segment_list;
};

struct cco_network* new_cco_network (struct user_options *options);

void build_segment (struct cco_network *the_network, const uint32_t index, const double new_pos[]);

void rescale_root (struct segment_node *iroot, const double Q_perf, const double delta_p);
void rescale_tree (struct segment_node *ibiff, const double Q_perf, const double delta_p, const int num_terminals);
void rescale_until_root (struct segment_node *ipar, struct segment_node *ipar_left, struct segment_node *ipar_right,\
                        const double Q_perf, const double delta_p, const int num_terminals);

void calc_middle_point_segment (struct segment_node *s, double pos[]);
void calc_relative_resistance_term (struct segment_node *iterm);
void calc_relative_resistance_subtree (struct segment_node *ibiff, struct segment_node *iconn, struct segment_node *inew);
void calc_radius_term (struct segment_node *iterm, const double Q_term, const double delta_p);
double calc_bifurcation_ratio (const double r1, const double r2, bool sign);

void check_bifurcation_rule (struct cco_network *the_network);

void grow_tree (struct cco_network *the_network);
void write_to_vtk (struct cco_network *the_network);

void test1 (struct cco_network *the_network);
void test2 (struct cco_network *the_network); 


// TODO: Move all this to a "utils.h"
// Auxiliar functions
double euclidean_norm (const double x1, const double y1, const double z1,\
                    const double x2, const double y2, const double z2);

#endif