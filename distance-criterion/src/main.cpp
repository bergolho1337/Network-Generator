// Author: Lucas Berg

#include <iostream>
#include <vector>
#include <algorithm>

#include <cstdio>
#include <cstdlib>
#include <cmath>

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

bool check_collision (Segment s1, Segment s2)
{
    Point *p1 = s1.p1;
    Point *p2 = s1.p2;
    Point *p3 = s2.p1;
    Point *p4 = s2.p2;

    double denominator = ((p2->x - p1->x) * (p3->y - p4->y)) - ((p2->y - p1->y) * (p3->x - p4->x));
    double numerator1 = ((p1->y - p4->y) * (p3->x - p4->x)) - ((p1->x - p4->x) * (p3->y - p4->y));
    double numerator2 = ((p1->y - p4->y) * (p2->x - p1->x)) - ((p1->x - p4->x) * (p2->y - p1->y));
    
    // Detect coincident lines (has a problem, read below)
    if (denominator == 0) return numerator1 == 0 && numerator2 == 0;

    double r = numerator1 / denominator;
    double s = numerator2 / denominator;

    bool intersect = (r >= 0 && r <= 1) && (s >= 0 && s <= 1);

    return intersect;
}

double euclidean_norm (const double x1, const double y1, const double z1,\
                       const double x2, const double y2, const double z2)
{
    return sqrt( pow(x2-x1,2) + pow(y2-y1,2) + pow(z2-z1,2) );
}

double calc_projection (Segment segment, Point inew)
{
    Point *distal = segment.p1;
    Point *prox = segment.p2;

    double length = euclidean_norm(prox->x,prox->y,prox->z,\
                                   distal->x,distal->y,distal->z);

    double dot_product = (prox->x - distal->x)*(inew.x - distal->x) +\
                         (prox->y - distal->y)*(inew.y - distal->y);
    
    return dot_product * pow(length,-2.0);
}

double calc_critical (Segment segment, Point inew, const double d_proj)
{
    Point *distal = segment.p1;
    Point *prox = segment.p2;

    if (d_proj >= 0.0 && d_proj <= 1.0)
    {
        double length = euclidean_norm(prox->x,prox->y,prox->z,\
                                   distal->x,distal->y,distal->z);

        double dot_product = (-prox->y + distal->y)*(inew.x - distal->x) +\
                         (prox->x - distal->x)*(inew.y - distal->y);

        return fabs(dot_product) * pow(length,-1.0);
    }
    else
    {
        double d_dist = euclidean_norm(inew.x,inew.y,inew.z,\
                                       distal->x,distal->y,distal->z);

        double d_prox = euclidean_norm(inew.x,inew.y,inew.z,\
                                       prox->x,prox->y,prox->z);
        
        return min(d_dist,d_prox);
    }
     
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
    segments.push_back(s1);
    print_segments(segments);

    double d_proj1 = calc_projection(segments[0],points[2]);
    double d_crit1 = calc_critical(segments[0],points[2],d_proj1);

    double d_proj2 = calc_projection(segments[0],points[3]);
    double d_crit2 = calc_critical(segments[0],points[3],d_proj2);

    double d_proj3 = calc_projection(segments[0],points[4]);
    double d_crit3 = calc_critical(segments[0],points[4],d_proj3);

    double d_proj4 = calc_projection(segments[0],points[5]);
    double d_crit4 = calc_critical(segments[0],points[5],d_proj4);

    printf("d_proj1 = %g\n",d_proj1);
    printf("d_crit1 = %g\n",d_crit1);

    printf("d_proj2 = %g\n",d_proj2);
    printf("d_crit2 = %g\n",d_crit2);
    
    printf("d_proj3 = %g\n",d_proj3);
    printf("d_crit3 = %g\n",d_crit3);

    printf("d_proj4 = %g\n",d_proj4);
    printf("d_crit4 = %g\n",d_crit4);

    write_to_vtk(points,segments);

    return 0;
}