//
// Created by bergolho on 12/02/19.
//

#ifndef CCO_H
#define CCO_H

#include <cstdio>
#include <cstdlib>
#include <cstdbool>
#include <cstdint>

#include "../point-list/point-list.h"
#include "../segment-list/segment-list.h"
#include "../options/user_options.h"

struct cco_network
{
    int num_terminals;

    int N_term;
    double Q_perf;
    double p_perf;
    double r_perf;
    double r_supp;

    struct point_list *point_list;
    struct segment_list *segment_list;
};

struct cco_network* new_cco_network (struct user_options *options);

void build_segment (struct cco_network *the_network, const uint32_t index, const double new_pos[]);

void calc_middle_point_segment (struct segment_node *s, double pos[]);

void grow_tree (struct cco_network *the_network);
void write_to_vtk (struct cco_network *the_network);

void test1 (struct cco_network *the_network);
void test2 (struct cco_network *the_network); 


#endif