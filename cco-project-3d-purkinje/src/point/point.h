//
// Created by bergolho on 12/02/19.
//

#ifndef POINT_H
#define POINT_H

#include <cstdio>
#include <cstdlib>
#include <cstdbool>
#include <cstdint>

class Point
{
public:
    uint32_t id;
    double x, y, z;
    double lat;
    bool is_active;
public:
    Point () { };
    Point (const uint32_t id) { this->id = id; }
    Point (const uint32_t id, const double pos[]) { this->id = id; this->x = pos[0]; this->y = pos[1]; this->z = pos[2]; this->lat = 0.0; this->is_active = false; }
    Point (Point *input) { this->id = input->id; this->x = input->x; this->y = input->y; this->z = input->z; this->lat = input->lat; this->is_active = input->is_active; }
    void setId (const uint32_t id) { this->id = id; }
    void setCoordinate (const double pos[]) { this->x = pos[0]; this->y = pos[1]; this->z = pos[2]; }
    void setLAT (const double lat) { this->lat = lat; }
    void setActive (const bool is_active) { this->is_active = is_active; }
    void print () { printf("id (x y z) [LAT] {Active} = %u (%g %g %g) [%g] {%d}\n",this->id,this->x,this->y,this->z,this->lat,(int)this->is_active); }
};

#endif