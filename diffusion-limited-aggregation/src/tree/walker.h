//
// Created by bergolho on 26/01/20.
//

#ifndef WALKER_H
#define WALKER_H

#include <cstdio>
#include <cstdlib>
#include <cstdbool>
#include <cstdint>
#include <cmath>
#include <vector>
#include <string>

#include "../utils/utils.h"
#include "../options/user_options.h"

const double RADIUS = 6.0; // Radius of the walker

struct walker
{
    double pos[3];
    bool stuck;
    double radius;
};

struct walker_node
{
    uint32_t id;
    struct walker *value;           
    struct walker_node *next;   
};

struct walker_list
{
    uint32_t num_nodes;
    struct walker_node *list_nodes;
};

struct walker* new_walker (struct user_options *the_options);
struct walker* new_walker (const double x, const double y, const double z);
void free_walker (struct walker *the_walker);

uint32_t is_stuck (struct walker_list *the_tree, struct walker *the_other);

void print_walker (struct walker *the_walker);



struct walker_list* new_walker_list ();
void free_walker_list (struct walker_list *l);

struct walker_node* insert_walker_node (struct walker_list *l, struct walker *s);
void delete_node (struct walker_list *l, const uint32_t index);
struct walker_node* search_walker_node (struct walker_list *l, const uint32_t index);
bool is_empty (struct walker_list *l);
void print_list (struct walker_list *l);
void write_list (struct walker_list *l, const uint32_t iter);

void order_list (struct walker_list *l);

struct walker_node* new_walker_node (uint32_t id, struct walker *s);
void free_walker_node (struct walker_node *n);


#endif