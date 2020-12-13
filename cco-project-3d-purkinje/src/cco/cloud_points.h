//
// Created by bergolho on 12/02/19.
//

#ifndef CLOUD_POINTS_H
#define CLOUD_POINTS_H

#include <cstdio>
#include <cstdlib>
#include <cstdbool>
#include <cstdint>
#include <cmath>

#include <vector>
#include <string>
#include <algorithm>

#include "../point/point.h"
#include "../utils/utils.h"

class Cloud_Point
{
public:
    std::string filename;
    std::vector<Point*> points;
public:
    Cloud_Point (std::string cloud_points_filename);
    ~Cloud_Point ();
};

#endif