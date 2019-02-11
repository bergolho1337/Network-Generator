//
// Created by bergolho on 11/02/19.
//

#ifndef LINKEDLIST_H
#define LINKEDLIST_H

#include <cstdio>
#include <cstdlib>
#include <cstdbool>
#include <cstdint>

struct node
{
    uint32_t id;
    void *value;
    struct node *next;   
};

struct linked_list
{
    uint32_t num_nodes;
    struct node *list_nodes;
};

#endif