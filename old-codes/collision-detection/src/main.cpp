// Author: Lucas Berg

#include <iostream>
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
    Segment s2(&points[2],&points[3]);
    segments.push_back(s2);
    print_segments(segments);

    bool has_collision = check_collision(segments[0],segments[1]);
    if (has_collision)
        printf("[!] Collision between segments!\n");
    else
        printf("[!] No collision!\n");

    write_to_vtk(points,segments);

    return 0;
}