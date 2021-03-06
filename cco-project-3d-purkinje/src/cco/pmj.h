#ifndef PMJ_DATA_H
#define PMJ_DATA_H

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>

#include <vector>
#include <string>
#include <list>
#include <algorithm>

#include "../point/point.h"
#include "../utils/utils.h"
#include "../options/pmj_config.h"
#include "../cost_function_library/custom_function.h"

//class Segment;
//class Point;

class PMJ
{
public:
    uint32_t package_size;
    uint32_t total_num_connected;
    uint32_t max_connection_tries;
    uint32_t connection_rate;
    uint32_t package_head;
    uint32_t cur_package;
    double region_radius;
    double lat_error_tolerance;
    std::string location_filename;
    std::vector<bool> connected;
    std::vector<double> error;
    std::vector<double> reference;
    std::vector<double> aprox;
    std::vector<Point*> points;
    std::list<uint32_t> package;
    std::vector<uint32_t> penalty;
    std::vector< std::vector<Point*> > packages;
    CostFunction *cost_fn;
public:
    PMJ ();
    PMJ (PMJConfig *config);
    ~PMJ ();
    void concatenate (PMJ *input);
    PMJ* copy ();
    void print ();
private:
    void initialize_parameters ();
    void initialize_parameters (PMJConfig *config);
    void initialize_arrays (PMJConfig *config);
    void initialize_points (PMJConfig *config);
    void initialize_packages ();
    void initialize_packages_2 ();
};

bool comparePoint (Point *a, Point *b);

#endif