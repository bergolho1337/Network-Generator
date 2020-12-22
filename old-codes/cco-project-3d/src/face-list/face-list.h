//
// Created by bergolho on 12/02/19.
//

#ifndef FACELIST_H
#define FACELIST_H

#include <cstdio>
#include <cstdlib>
#include <cstdbool>
#include <cstdint>

struct face
{
    double x1, y1, z1;
    double x2, y2, z2;
    double x3, y3, z3;
    double nx, ny, nz;
};

struct face* new_face (const double pos1[], const double pos2[], const double pos3[], const double n[]);
void free_face (struct face *f);

struct face_node
{
    uint32_t id;
    struct face *value;           
    struct face_node *next;   
};

struct face_list
{
    uint32_t num_nodes;
    struct face_node *list_nodes;
};

struct face_list* new_face_list ();
void free_face_list (struct face_list *l);

struct face_node* insert_point (struct face_list *l,\
                            const double pos1[], const double pos2[], const double pos3[], const double n[]);
void delete_node (struct face_list *l, const uint32_t index);
struct face_node* search_node (struct face_list *l, const uint32_t index);
bool is_empty (struct face_list *l);
void print_list (struct face_list *l);
void write_list (struct face_list *l, FILE *log_file); 

void order_list (struct face_list *l);

struct face_node* new_node (uint32_t id, const double pos1[], const double pos2[], const double pos3[], const double n[]);
void free_node (struct node *n);

#endif