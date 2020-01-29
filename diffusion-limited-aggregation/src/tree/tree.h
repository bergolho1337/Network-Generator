//
// Created by bergolho on 26/01/20.
//

#ifndef TREE_H
#define TREE_H

#include <cstdio>
#include <cstdlib>
#include <cstdbool>
#include <cstdint>
#include <vector>

#include "walker.h"
#include "segment.h"

#include "../options/user_options.h"

struct dla_tree
{
    struct walker_list *point_list;
    struct segment_list *segment_list;
};

struct dla_tree* new_dla_tree ();
void free_dla_tree (struct dla_tree *the_tree);

void grow_tree (struct dla_tree *the_tree, struct user_options *the_options);
void write_to_vtk (struct dla_tree *the_tree);

#endif