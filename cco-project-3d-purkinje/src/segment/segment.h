//
// Created by bergolho on 12/02/19.
//

#ifndef SEGMENT_H
#define SEGMENT_H

#include <cstdio>
#include <cstdlib>
#include <cstdbool>
#include <cstdint>
#include <cmath>

#include "../point/point.h"

class Segment
{
public:
    uint32_t id;
    double radius;
    double beta;
    double length;
    //double resistance;
    //double Q;
    //double delta_p;
    int ndist;

    bool can_touch;
    bool prune;

    Point *src;
    Point *dest;

    Segment *left;
    Segment *right;
    Segment *parent;
public:
    Segment (const uint32_t id, Point *src, Point *dest, Segment *left, Segment *right, Segment *parent)
    {
        this->id = id;
        this->src = src;
        this->length = sqrt( powf(src->x-dest->x,2) + powf(src->y-dest->y,2) + powf(src->z-dest->z,2) );
        this->dest = dest;
        this->left = left;
        this->right = right;
        this->parent = parent;
        this->ndist = 1;
        this->beta = -1;
        this->radius = 0.0;
        this->prune = false;
        this->can_touch = true;
    }
    ~Segment ()
    {
        this->parent = NULL;
        this->left = NULL;
        this->right = NULL;
        this->src = NULL;
        this->dest = NULL;
    }
    void calc_middle_point (double pos[])
    {
        Point *src = this->src;
        Point *dest = this->dest;

        pos[0] = (src->x + dest->x) / 2.0;
        pos[1] = (src->y + dest->y) / 2.0;
        pos[2] = (src->z + dest->z) / 2.0;        
    }
    void calc_unitary_vector (double u[])
    {
        if (this->length < 1.0e-08) this->length = 1.0e-08;

        u[0] = (this->dest->x - this->src->x) / this->length;
        u[1] = (this->dest->y - this->src->y) / this->length;
        u[2] = (this->dest->z - this->src->z) / this->length;
    }
    double calc_radius ()
    {
        if (this->parent == NULL)
            return this->radius;
        else
            return this->beta * this->parent->calc_radius();
    }
    double calc_length ()
    {
        return sqrt( powf(src->x-dest->x,2) + powf(src->y-dest->y,2) + powf(src->z-dest->z,2) );
    }
    double calc_middle_point_lat ()
    {
        return (src->lat + dest->lat)/2.0;
    }
    uint32_t calc_level ()
    {
        uint32_t result = 0;
        Segment *tmp = this;
        while (tmp != NULL)
        {
            result++;
            tmp = tmp->parent;
        }
        return result;
    }
    bool is_terminal ()
    {
        return (this->left == NULL && this->right == NULL) ? true : false;
    }
    bool is_inside_region (const double center[], const double radius)
    {
        double dist = sqrt( powf(dest->x-center[0],2) + powf(dest->y-center[1],2) + powf(dest->z-center[2],2) ); 
        return (dist < radius) ? true : false;
    }
    bool is_touchable ()
    {
        return this->can_touch;
    }
    void print ()
    {
        printf("Segment %u (%u,%u) (%g %g %g) - (%g %g %g) -- NDIST = %d -- RADIUS = %g mm\n",this->id,src->id,dest->id,\
                                                                            src->x,src->y,src->z,dest->x,dest->y,dest->z,\
                                                                            this->ndist,this->radius);
        printf("\tBETA = %g\n",this->beta);

        if (this->parent == NULL)
            printf("\tPARENT = NIL\n");
        else
            printf("\tPARENT = %u\n",this->parent->id);
        if (this->left == NULL)
            printf("\tLEFT = NIL\n");
        else
            printf("\tLEFT = %u\n",this->left->id);
        if (this->right == NULL)
            printf("\tRIGHT = NIL\n");
        else
            printf("\tRIGHT = %u\n",this->right->id);
    }
};

#endif