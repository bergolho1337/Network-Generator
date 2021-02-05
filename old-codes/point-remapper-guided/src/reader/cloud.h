#ifndef CLOUD_H_
#define CLOUD_H_

#include <iostream>
#include <vector>
#include <cstdint>

using namespace std;

class Point_3d
{
public:
    uint32_t id;
    double x, y, z;
public:
    Point_3d (const uint32_t id, const double x, const double y, const double z)
    {
        this->id = id;
        this->x = x;
        this->y = y;
        this->z = z;
    }
};

class Cloud_Point
{
public:
    vector<bool> taken;
    vector<Point_3d> the_points;
    vector<Point_3d> the_remapped_points;
public:
    Cloud_Point (const char filename[]);
    void print ();
};

#endif
