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
    std::vector<double> error;
    std::vector<double> ref;
    std::vector<double> aprox;
public:
    Cloud ();
    Cloud (CloudConfig *config);
    ~Cloud ();
    uint32_t sort_point (Point *p, const uint32_t rand_offset);
    Cloud* copy ();
    Cloud* get_points_around_region (Point *pmj_point, const uint32_t num_points, const double region_radius);
    void concatenate (Cloud *input);
    void print ();
};

#endif