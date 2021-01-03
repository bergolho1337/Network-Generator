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
#include "../cco/constants.h"

class Segment
{
public:
    uint32_t id;
    double radius;
    double beta;
    double length;
    int ndist;

    bool can_touch;
    bool prune;
    bool adjusted;

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
        this->adjusted = false;
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
    double calc_propagation_velocity ()
    {
        const double G = 7.9;
        const double Cf = 3.4;
        const double tauf = 0.1;

        double d = this->radius*2.0 * 1000.0;   // {mm}->{um}
    
        // Output in {m/s}
        return powf( (G*d)/(4.0*Cf*tauf) , 0.5 ) * 0.1;
    }
    double calc_pathway_length ()
    {
        double result = 0.0;
        Segment *tmp = this;
        while (tmp != NULL)
        {
            result += tmp->length;
            tmp = tmp->parent;
        }
        return result;
    }
    double calc_terminal_local_activation_time ()
    {
        double result = 0.0;
        Segment *tmp = this;
        while (tmp != NULL)
        {
            double cv = tmp->calc_propagation_velocity()*M_S_TO_UM_MS;    // {m/s}->{um/ms}
            double dist = tmp->length;
            double lat = dist*M_TO_UM/cv;

            result += lat;
            
            tmp = tmp->parent;
        }

        // Return the LAT in {ms}
        return result; 
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
    void update_radius (const double cv)
    {
        const double G = 7.9;
        const double Cf = 3.4;
        const double tauf = 0.1;

        double cv_m_per_s = cv / 1000.0;
        double new_diameter = (cv_m_per_s*cv_m_per_s*4.0*Cf*tauf) / (G) * 100.0; // {um}
        new_diameter *= 1.0e-03;

        this->radius = new_diameter/2.0;
        this->adjusted = true;
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

void eliminate_segment_from_list (std::vector<Segment*> &s_list, Segment *s);
void order_segment_list (std::vector<Segment*> &s_list);

#endif