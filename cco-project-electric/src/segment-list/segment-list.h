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
    double I;             // Current
    double delta_v;       // Potential drop 
    double radius;        // Radius 
    double beta;          // Relative radius (ratio between my radius and parent radius = beta)
    double resistance;    // Relative resistance R*
    double length;
    uint32_t ndist;

    struct point_node *src;
    struct point_node *dest;

    struct segment_node *left;      
    struct segment_node *right;     
    struct segment_node *parent;

};

struct segment* new_segment (struct point_node *src, struct point_node *dest,\
                        struct segment_node *left, struct segment_node *right, struct segment_node *parent,\
                        const double I, const double v);
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

struct segment_node* insert_segment_node (struct segment_list *l, struct segment *s);
void delete_node (struct segment_list *l, const uint32_t index);
struct segment_node* search_segment_node (struct segment_list *l, const uint32_t index);
bool is_empty (struct segment_list *l);
void print_list (struct segment_list *l);
void write_list (struct segment_list *l, FILE *log_file); 

void order_list (struct segment_list *l);

struct segment_node* new_segment_node (uint32_t id, struct segment *s);
void free_segment_node (struct segment_node *n);

#endif