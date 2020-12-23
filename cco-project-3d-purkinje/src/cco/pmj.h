#ifndef PMJ_DATA_H
#define PMJ_DATA_H

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>

#include <vector>
#include <string>
#include <algorithm>

#include "../point/point.h"
#include "../utils/utils.h"
#include "../options/pmj_config.h"
#include "../cost_function_library/cost_function.h"

//class Segment;
//class Point;

class PMJ
{
public:
    uint32_t total_num_connected;
    uint32_t max_connection_tries;
    uint32_t connection_rate;
    double region_radius;
    double lat_error_tolerance;
    std::string location_filename;
    std::vector<bool> connected;
    std::vector<double> error;
    std::vector<double> aprox;
    std::vector<Point*> points;
    CostFunction *cost_fn;
public:
    PMJ ();
    PMJ (PMJConfig *config);
    ~PMJ ();
    void print ();
    PMJ* copy ();
};

bool comparePoint (Point *a, Point *b);

#endif