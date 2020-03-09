// Author: Lucas Berg
// Program that reads a STL file and store its faces and point on a vector.

#ifndef READER_H
#define READER_H

#include <iostream>
#include <cmath>
#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <map>

using namespace std;

class Point_Table
{
public:
    double x, y, z;
public:
    Point_Table () {};
    Point_Table (const double x, const double y, const double z)
    {
        this->x = x;
        this->y = y;
        this->z = z;
    };
    void print ()
    {
        printf("(%.2lf,%.2lf,%.2lf)\n",this->x,this->y,this->z);
    };
    bool operator <(const Point_Table& pt) const
    {
        return (this->x < pt.x) || \
               ((!(pt.x < this->x)) && (this->y < pt.y)) || \
               ((!(pt.x < this->x)) && (!(pt.y < this->y)) && (this->z < pt.z));
    }
};

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
    bool operator <(const Point& pt) const
    {
        return (this->x < pt.x) || \
               ((!(pt.x < this->x)) && (this->y < pt.y)) || \
               ((!(pt.x < this->x)) && (!(pt.y < this->y)) && (this->z < pt.z));
    }
};

class Face
{
public:
    Point *v1;
    Point *v2;
    Point *v3;
    double normal[3];
public:
    Face (Point *a, Point *b, Point *c, double n[])
    {
        this->v1 = a;
        this->v2 = b;
        this->v3 = c;
        for (uint32_t i = 0; i < 3; i++)
            this->normal[i] = n[i];
    };
    void print ()
    {
        printf("\tVertex = [%lf %lf %lf] [%lf %lf %lf] [%lf %lf %lf] -- {%u %u %u}\n",this->v1->x,this->v1->y,this->v1->z,\
                                                        this->v2->x,this->v2->y,this->v2->z,\
                                                        this->v3->x,this->v3->y,this->v3->z,\
                                                        this->v1->id,this->v2->id,this->v3->id);
        printf("\tNormal = (%lf %lf %lf)\n",this->normal[0],this->normal[1],this->normal[2]);
    };
};

void read_face (FILE *file, vector<Face> &faces);
void read_faces (const char *filename, vector<Face> &faces, map<Point_Table,uint32_t> &point_table);
void insert_point_into_table (map<Point_Table,uint32_t> &point_table, const double pos[]);
void print_points (vector<Point> points);
void print_faces (vector<Face> faces);

#endif