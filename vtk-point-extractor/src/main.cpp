// Author: Lucas Berg

#include <iostream>
#include <cmath>
#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <vector>

using namespace std;

class Point
{
public:
    int id;
    double x, y, z;
public:
    Point () {};
    Point (const int id, const double x, const double y, const double z)
    {
        this->id = id;
        this->x = x;
        this->y = y;
        this->z = z;
    };
    void print ()
    {
        printf("%d -- (%.2lf,%.2lf,%.2lf)\n",this->id,this->x,this->y,this->z);
    };
};

void read_points (const char *filename, vector<Point> &points)
{
    char str[200];
    uint32_t num_points;
    double x, y, z;

    FILE *file = fopen(filename,"r");

    while (fscanf(file,"%s",str) != EOF)
    {
        if (strcmp(str,"POINTS") == 0)
            break;
    }
    fscanf(file,"%u %s",&num_points,str);
    
    for (uint32_t i = 0; i < num_points; i++)
    {
        fscanf(file,"%lf %lf %lf",&x,&y,&z);

        Point point(i,x,y,z);
        points.push_back(point);
    }

    fclose(file);
}

void rescale_points (vector<Point> &points, const double scale_factor)
{
    for (uint32_t i = 0; i < points.size(); i++)
    {
        points[i].x *= scale_factor;
        points[i].y *= scale_factor;
        points[i].z *= scale_factor;
    }
}

void print_points (vector<Point> points)
{
    for (int i = 0; i < (int)points.size(); i++)
        points[i].print();
}

void write_points_to_pts (vector<Point> points)
{
    FILE *file = fopen("outputs/points.pts","w+");

    fprintf(file,"%lu\n",points.size());
    for (uint32_t i = 0; i < points.size(); i++)
        fprintf(file,"%lf %lf %lf\n",points[i].x,points[i].y,points[i].z);
    
    fclose(file);
}

int main (int argc, char *argv[])
{
    if (argc-1 != 1)
    {
        printf("Usage:> %s <input_filename>\n",argv[0]);
        exit(EXIT_FAILURE);   
    }

    vector<Point> points;
    read_points(argv[1],points);
    //print_points(points);

    double scale_factor = 0.025;
    rescale_points(points,scale_factor);

    write_points_to_pts(points);

    return 0;
}
