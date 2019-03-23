#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <vector>

#define PRINT_LINE "================================================================="
#define PRINT_LINE_2 "------------------------------------------------------------------"

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