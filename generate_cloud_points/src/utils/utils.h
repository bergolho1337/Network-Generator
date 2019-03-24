#ifndef UTILS_H
#define UTILS_H

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <map>

// *********************************************************************************************
// CONSTANTS AND MACROS
#define PRINT_LINE "================================================================="
#define PRINT_LINE_2 "------------------------------------------------------------------"

// *********************************************************************************************

class Point
{
public:
    double x, y, z;
public:
    Point () {};
    Point (const double x, const double y, const double z)
    {
        this->x = x;
        this->y = y;
        this->z = z;
    };
    void print ()
    {
        printf("(%.2lf,%.2lf,%.2lf)\n",this->x,this->y,this->z);
    };
};

void usage (const char pname[]);
double generate_random_number (); 
double calc_euclidean_dist (const double a[], const double b[]);
bool check_point (const double new_pos[], std::vector<Point> points, const double tolerance);
bool get_parameter_value_from_map (std::map<std::string,double> *params,\
                                const std::string key, double *value);

#endif