#ifndef FRACTAL_H
#define FRACTAL_H

#include <cstdio>
#include <cstdlib>
#include <cstdbool>
#include <cstdint>
#include <cmath>

#include <vector>
#include <queue>

#include "../utils/utils.h"

class Point
{
public:
    uint32_t id;
    double x, y, z;
public:
    Point (const uint32_t id, const double x, const double y, const double z)
    {
        this->id = id;
        this->x = x;
        this->y = y;
        this->z = z;
    }
    void print ()
    {
        printf("Point %u --> (%g,%g,%g)\n",this->id,this->x,this->y,this->z);
    }
};

class Line
{
public:
    uint32_t id;
    uint32_t src;
    uint32_t dest;
    double diameter;
public:
    Line (const uint32_t id, const uint32_t src_id, const uint32_t dest_id, const double diameter)
    {
        this->id = id;
        this->src = src_id;
        this->dest = dest_id;
        this->diameter = diameter;
    }
    void print ()
    {
        printf("Line %u --> || src = %u || dest = %u || diameter = %g ||\n",this->id,this->src,this->dest,this->diameter);
    }
};

class Fractal_Tree
{
public:
    uint32_t max_iterations;
    double root_pos[3];

    double initial_length;
    double initial_angle;
    double initial_diameter;

    double angle_decrease_ratio;
    double length_decrease_ratio;
    double diameter_decrease_ratio;

    std::vector<Point> the_points;
    std::vector<Line> the_lines; 
public:
    Fractal_Tree ();
    void grow_network ();
    void make_root (std::queue<Line> &q);
    void print ();
    void write ();
};


#endif
