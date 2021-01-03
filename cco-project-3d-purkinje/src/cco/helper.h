#ifndef HELPER_H
#define HELPER_H

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>

#include <vector>
#include <string>
#include <list>
#include <algorithm>

#include "../segment/segment.h"
#include "../utils/utils.h"
#include "../options/local_optimization_config.h"

#include "constants.h"

bool check_null_segments (std::vector<Segment*> s_list);
bool check_bifurcation_rule (FILE *log_file, const double gamma, std::vector<Segment*> s_list);
bool has_collision (std::vector<Segment*> s_list, Segment *s, Point *p);
bool distance_criterion (Segment *s, const double pos[], const double d_threash);
bool connection_search (std::vector<Segment*> s_list, Point *p, const double d_threash);
bool check_collisions_and_fill_feasible_segments (std::vector<Segment*> s_list, Point *p, std::vector<Segment*> &feasible_segments);
void update_segment_pointers (Segment *iconn, Segment *ibiff, Segment *inew, const bool is_restore);
void get_segment_length (std::vector<Segment*> s_list, std::vector<double> &segments);
void get_bifurcation_angles(std::vector<Segment*> s_list, std::vector<double> &angles);
Segment* get_terminal (std::vector<Point*> p_list, std::vector<Segment*> s_list, const double pos[]);
uint32_t update_ndist (Segment *s, const bool is_restore);
Point* generate_bifurcation_node (std::vector<Point*> p_list, Segment *iconn, LocalOptimizationConfig *local_opt_config, const bool using_local_optimization);
Point* generate_terminal_node (std::vector<Point*> p_list, Point *p);

#endif