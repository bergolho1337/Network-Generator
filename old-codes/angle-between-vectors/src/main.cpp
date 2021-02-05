// Author: Lucas Berg

#include <iostream>
#include <vector>
#include <algorithm>

#include <cstdio>
#include <cstdlib>
#include <cmath>

using namespace std;

double euclidean_norm (const double x1, const double y1, const double z1,\
                       const double x2, const double y2, const double z2);

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

class Segment
{
public:
    Point *p1;
    Point *p2;
public:
    Segment (Point *a, Point *b)
    {
        this->p1 = a;
        this->p2 = b;
    };
    void print ()
    {
        printf("(%d -- %d)\n",this->p1->id,this->p2->id);
    };
    double norm ()
    {
        return euclidean_norm(this->p1->x,this->p1->y,this->p1->z,\
                              this->p2->x,this->p2->y,this->p2->z);
    }
};

void read_points (const char *filename, vector<Point> &points)
{
    FILE *file = fopen(filename,"r");
    double x, y, z;

    while (fscanf(file,"%lf %lf %lf",&x,&y,&z) != EOF)
    {
        Point p(points.size(),x,y,z);
        points.push_back(p);
    }
    fclose(file);
}

void print_points (vector<Point> points)
{
    for (int i = 0; i < (int)points.size(); i++)
        points[i].print();
}

void print_segments (vector<Segment> segments)
{
    for (int i = 0; i < (int)segments.size(); i++)
        segments[i].print();
}

void write_to_vtk (vector<Point> points, vector<Segment> segments)
{
    FILE *file = fopen("lines.vtk","w+");

    fprintf(file,"# vtk DataFile Version 3.0\n");
    fprintf(file,"Line Collision\n");
    fprintf(file,"ASCII\n");
    fprintf(file,"DATASET POLYDATA\n");
    fprintf(file,"POINTS %lu float\n",points.size());
    for (unsigned int i = 0; i < points.size(); i++)
        fprintf(file,"%g %g %g\n",points[i].x,points[i].y,points[i].z);
    fprintf(file,"LINES %lu %lu\n",segments.size(),segments.size()*3);
    for (unsigned int i = 0; i < segments.size(); i++)
        fprintf(file,"2 %d %d\n",segments[i].p1->id,segments[i].p2->id);

    fclose(file);
}

double euclidean_norm (const double x1, const double y1, const double z1,\
                       const double x2, const double y2, const double z2)
{
    return sqrt( pow(x2-x1,2) + pow(y2-y1,2) + pow(z2-z1,2) );
}

double dot_product (Segment u, Segment v)
{
    return (u.p2->x - u.p1->x)*(v.p2->x - v.p1->x) +\
           (u.p2->y - u.p1->y)*(v.p2->y - v.p1->y) +\
           (u.p2->z - u.p1->z)*(v.p2->z - v.p1->z); 
}

double calc_angle_between_vectors (Segment u, Segment v)
{
    double norm_u = u.norm();
    double norm_v = v.norm();

    double uv = dot_product(u,v);
    
    return acos(uv / (norm_u * norm_v));
}

double convert_radians_to_degree (const double radians)
{
    return radians * 180.0 / M_PI;
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
    print_points(points);

    vector<Segment> segments;
    Segment s1(&points[0],&points[1]);
    Segment s2(&points[2],&points[3]);
    segments.push_back(s1);
    segments.push_back(s2);
    print_segments(segments);

    double angle_radians = calc_angle_between_vectors(segments[0],segments[1]);
    printf("Angle = %g degrees\n",convert_radians_to_degree(angle_radians));

    write_to_vtk(points,segments);

    return 0;
}