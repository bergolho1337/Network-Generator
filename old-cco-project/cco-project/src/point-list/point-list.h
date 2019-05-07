//
// Created by bergolho on 12/02/19.
//

#ifndef POINTLIST_H
#define POINTLIST_H

#include <cstdio>
#include <cstdlib>
#include <cstdbool>
#include <cstdint>

struct point
{
    double x, y, z;
};

struct point* new_point (const double pos[]);
void free_point (struct point *p);

struct point_node
{
    uint32_t id;
    struct point *value;           
    struct point_node *next;   
};

struct point_list
{
    uint32_t num_nodes;
    struct point_node *list_nodes;
};

struct point_list* new_point_list ();
void free_point_list (struct point_list *l);

struct point_node* insert_point (struct point_list *l, const double pos[]);
void delete_node (struct point_list *l, const uint32_t index);
struct point_node* search_node (struct point_list *l, const uint32_t index);
bool is_empty (struct point_list *l);
void print_list (struct point_list *l);
void write_list (struct point_list *l, FILE *log_file); 

void order_list (struct point_list *l);

struct point_node* new_node (uint32_t id, const double pos[]);
void free_node (struct node *n);

#endif