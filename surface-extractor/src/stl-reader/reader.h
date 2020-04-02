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

#include <vtkTriangle.h>
#include <vtkPointData.h>
#include <vtkSTLReader.h>
#include <vtkSTLWriter.h>
#include <vtkPolyDataMapper.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>

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
    uint32_t vertex_index_1;
    uint32_t vertex_index_2;
    uint32_t vertex_index_3;

    bool marked;
public:
    Face (const uint32_t id_1, const uint32_t id_2, const uint32_t id_3)
    {
        this->vertex_index_1 = id_1;
        this->vertex_index_2 = id_2;
        this->vertex_index_3 = id_3;

        this->marked = false;
    };
    void print ()
    {
        printf("\tVertexes = {%u %u %u}\n",this->vertex_index_1,this->vertex_index_2,this->vertex_index_3);
    };
};

void read_faces (const char *filename, vector<Point> &points, vector<Face> &faces);
void read_points_from_pts (const char *filename, vector<Point> &points);
void insert_point_into_table (map<Point_Table,uint32_t> &point_table, const double pos[]);
void print_points (vector<Point> points);
void print_faces (vector<Face> faces);
void write_marked_faces_to_stl(vector<Point> points, vector<Face> faces);
void write_face_to_stl(FILE *file, Face face, vector<Point> points);

#endif