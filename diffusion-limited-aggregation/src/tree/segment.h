//
// Created by bergolho on 26/01/20.
//

#ifndef SEGMENT_H
#define SEGMENT_H

#include <cstdio>
#include <cstdlib>
#include <cstdbool>
#include <cstdint>
#include <vector>

struct segment
{
    uint32_t src;
    uint32_t dest;
};

struct segment* new_segment (const uint32_t src, const uint32_t dest);
void free_segment (struct segment *the_segment);

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