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
#include "../options/cloud_config.h"

class Cloud
{
public:
    std::string filename;
    uint32_t cur_index;
    std::vector<Point*> points;
    std::vector<bool> connected;
public:
    Cloud ();
    Cloud (CloudConfig *config);
    ~Cloud ();
    uint32_t sort_point (Point *p, const uint32_t rand_offset);
    Cloud* copy ();
    void concatenate (Cloud *input);
    void print ();
};

#endif