//
// Created by bergolho on 11/02/19.
//

#ifndef LINKEDLIST_H
#define LINKEDLIST_H

#include <cstdio>
#include <cstdlib>
#include <cstdbool>
#include <cstdint>

// Change the node datatype of list ... HERE !!!
#define NODE_DATATYPE double

struct node
{
    uint32_t id;
    NODE_DATATYPE value;           
    struct node *next;   
};

struct linked_list
{
    uint32_t num_nodes;
    struct node *list_nodes;
};

struct linked_list* new_linked_list ();
void free_linked_list (struct linked_list *l);

void insert_node (struct linked_list *l, const NODE_DATATYPE key);
void delete_node (struct linked_list *l, const NODE_DATATYPE key);
struct node* search_node (struct linked_list *l, const NODE_DATATYPE key);
bool is_empty (struct linked_list *l);
void print_list (struct linked_list *l); 

void order_list (struct linked_list *l);

struct node* new_node (uint32_t id, const NODE_DATATYPE key);
void free_node (struct node *n);

#endif