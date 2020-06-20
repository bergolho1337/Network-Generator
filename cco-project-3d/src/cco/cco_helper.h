//
// Created by bergolho on 12/02/19.
//

#ifndef CCO_HELPER_H
#define CCO_HELPER_H

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

void calc_middle_point_segment (struct segment_node *s, double pos[]);
void calc_unitary_vector (struct segment_node *s, double u[]);
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
                        const double G, const double Cf, const double tau_f);
double calc_total_activation_time (struct cco_network *the_network,\
                            const double G, const double Cf, const double tau_f);
double calc_segment_activation_time (struct segment_node *s,\
                        const double G, const double Cf, const double tau_f);
double calc_propagation_velocity (const double d,\
                        const double G, const double Cf, const double tau_f);
double calc_lambda_m (const double r, const double rc, const double rm);
double calc_tau_m (const double cm, const double rm);

double calc_segment_activation_time_using_level (const double at, struct segment_node *iconn);
double calc_segment_level (struct segment_node *iconn);

double calc_custom_function (struct cco_network *the_network, const double beta, const double alpha);
double calc_segment_custom_function (struct segment_node *s, const double beta, const double alpha);
double calc_segment_custom_function_with_level_penalty (const double eval, struct segment_node *iconn);


bool has_deviation (struct segment_list *s_list, struct segment_node *inew,\
                    const double new_at, const double limit,\
                    const double G, const double Cf, const double tau_f);
bool check_angle_restriction (const double angle, const double min_angle, const double max_angle);
bool is_terminal (struct segment_node *s);

#endif