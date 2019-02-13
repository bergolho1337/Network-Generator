//
// Created by bergolho on 12/02/19.
//

#ifndef SEGMENTLIST_H
#define SEGMENTLIST_H

#include <cstdio>
#include <cstdlib>
#include <cstdbool>
#include <cstdint>

#include "../point-list/point-list.h"

struct segment
{
    double beta_l;
    double beta_r;
    double Q;
    double p;
    double radius;
    double resistance;          // Relative resistance
    double length;
    uint32_t ndist;

    struct point_node *src;
    struct point_node *dest;

    struct segment *left;
    struct segment *right;
    struct segment *parent;

};

struct segment* new_segment (struct point_node *src, struct point_node *dest,\
                        struct segment *left, struct segment *right, struct segment *parent,\
                        const double Q, const double p);
void free_segment (struct segment *s);

struct segment_node
{
    uint32_t id;
    struct segment *value;           
    struct segment_node *next;   
};

struct segment_list
{
    uint32_t num_nodes;
    struct segment_node *list_nodes;
};

struct segment_list* new_segment_list ();
void free_segment_list (struct segment_list *l);

void insert_segment_node (struct segment_list *l, struct segment *s);
void delete_node (struct segment_list *l, const uint32_t index);
struct segment_node* search_node (struct segment_list *l, const uint32_t index);
bool is_empty (struct segment_list *l);
void print_list (struct segment_list *l); 

void order_list (struct segment_list *l);

struct segment_node* new_segment_node (uint32_t id, struct segment *s);
void free_segment_node (struct segment_node *n);

#endif