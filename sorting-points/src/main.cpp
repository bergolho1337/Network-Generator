// Author: Lucas Berg

#include <iostream>
#include <vector>
#include <random>
#include <algorithm>

#include <cstdio>
#include <cstdlib>
#include <cmath>

using namespace std;

// ====================================================================================================================
// CONSTANTS AND MACROS
#define PRINT_LINE "--------------------------------------------------------------------------------------------"

class Point
{
public:
    double x, y, z;
public:
    Point (const double x, const double y, const double z)
    {
        this->x = x;
        this->y = y;
        this->z = z;
    }
    void print ()
    {
        printf("(%g,%g,%g)\n",this->x,this->y,this->z);
    }
};

// Sort points first by 'x', then by 'y' and finally by 'z'
bool compare (const Point &a, const Point &b)
{
    if (a.x < b.x)
        return true;
    else
    {
        if (a.x == b.x)
        {
            if (a.y < b.y)
                return true;
            else
            {
                if (a.y == b.y)
                {
                    if (a.z < b.z)
                        return true;
                    else
                        return false;
                }
                else
                {
                    return false;
                }
            }
        }
        else
        {
            return false;
        }
    }
}

// ====================================================================================================================

void read_user_input (vector<Point> &points, const char filename[])
{
    FILE *file = fopen(filename,"r");

    uint32_t np;
    fscanf(file,"%u",&np);

    for (uint32_t i = 0; i < np; i++)
    {
        double pos[3];
        fscanf(file,"%lf %lf %lf",&pos[0],&pos[1],&pos[2]);
        Point p(pos[0],pos[1],pos[2]);
        points.push_back(p);
    }

    fclose(file);
}

void print_points (vector<Point> points)
{
    printf("Number of points = %u\n",points.size());
    for (uint32_t i = 0; i < points.size(); i++)
        points[i].print();
}

int main (int argc, char *argv[])
{
    if (argc-1 != 1)
    {
        printf("%s\n",PRINT_LINE);
        printf("Usage:> %s <input_file>\n",argv[0]);
        printf("%s\n",PRINT_LINE);
        printf("<input_file> = Input file with the points coordinates\n");
        exit(EXIT_FAILURE);   
    }

    vector<Point> main_points;
    read_user_input(main_points,argv[1]);
    print_points(main_points);

    sort(main_points.begin(),main_points.end(),compare);

    print_points(main_points);

    return 0;
}
