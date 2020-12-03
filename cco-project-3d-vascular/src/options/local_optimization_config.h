#ifndef LOCAL_OPTIMIZATION_CONFIG_H
#define LOCAL_OPTIMIZATION_CONFIG_H

#include <cstdio>
#include <cstdlib>
#include <cmath>

#include "../options/cost_function_config.h"
#include "../cco/cco.h"

struct local_optimization_config;

// Flag to allow a function export ...
#define EXPORT_FN

// Template for a local optimization function
#define SET_LOCAL_OPTIMIZATION_FUNCTION(name) EXPORT_FN void name(struct segment_node *iconn,\
                                        struct segment_node *ibiff,\
                                        struct segment_node *inew,\
                                        std::vector<struct point> &test_positions)
typedef SET_LOCAL_OPTIMIZATION_FUNCTION(set_local_optimization_function_fn);

struct local_optimization_config
{
    void *handle;

    char *function_name;
    char *library_name;
    //std::map<std::string,double> *params;    // Parameters of the local optimization function

    bool first_call;
    double best_pos[3];                                    // Best position for the bifurcation

    set_local_optimization_function_fn *function;          // Reference to the local optimization function

};

struct local_optimization_config* new_local_optimization_config ();
void free_local_optimization_config (struct local_optimization_config *config);

void set_local_optimization_function (struct local_optimization_config *config);

void print_local_optimization_function_config (struct local_optimization_config *config);

// Auxiliary functions
void save_original_bifurcation_position (struct segment_node *iconn, double ori_pos[]);
void initialize_best_position_as_middle_point(double best_pos[], const double ori_pos[]);
bool is_corner (const uint32_t i, const uint32_t j, const uint32_t NE);
void move_bifurcation_location (struct segment_node *iconn, struct segment_node *ibiff, struct segment_node *inew,\
                            const double pos[]);

#endif