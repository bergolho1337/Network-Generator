#ifndef LOCAL_OPTIMIZATION_CONFIG_H
#define LOCAL_OPTIMIZATION_CONFIG_H

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <string>

class Segment;
class Point;

class LocalOptimizationConfig
{
public:
    bool first_call;
    std::string function_name;
    std::string library_name;
    double best_pos[3];
public:
    LocalOptimizationConfig ();
    ~LocalOptimizationConfig ();
    void initialize_best_position_as_middle_point (double best_pos[], const double ori_pos[]);
    void save_original_bifurcation_position (Segment *s, double ori_pos[]);
    void update_bifurcation_position (Point *p);
    void print ();
};

// Auxiliary functions
/*
void initialize_best_position_as_middle_point(double best_pos[], const double ori_pos[]);
bool is_corner (const uint32_t i, const uint32_t j, const uint32_t NE);
bool is_vertex (const uint32_t i, const uint32_t j, const uint32_t NE);
void move_bifurcation_location (struct segment_node *iconn, struct segment_node *ibiff, struct segment_node *inew,\
                            const double pos[]);
*/

#endif