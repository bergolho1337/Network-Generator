#ifndef UTILS_H
#define UTILS_H

#include <cstdio>
#include <cstdlib>
#include <cstdbool>
#include <cstdint>
#include <cstring>
#include <cmath>
#include <vector>

class Point
{
public:
    Point (const uint32_t id, const double x, const double y, const double z)
    {
        this->id = id;
        this->x = x;
        this->y = y;
        this->z = z;
        this->lat = 0;
        this->is_active = 0;
    }
public:
    uint32_t id;
    double x, y, z;
    double lat;
    int is_active;
};

double calc_norm (const double x1, const double y1, const double z1,\
                const double x2, const double y2, const double z2);
uint32_t calc_closest_point (Point p, std::vector<Point> points);


#endif