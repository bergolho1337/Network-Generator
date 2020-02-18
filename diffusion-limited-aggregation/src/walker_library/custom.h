#ifndef SQUARE_LIB_H
#define SQUARE_LIB_H

#include <cstdio>
#include <cstdlib>
#include <cmath>

#include <map>
#include <vector>

#include "../options/walker_config.h"
#include "../tree/tree.h"

class Point_Custom
{
public:
    double x, y, z;
public:
    Point_Custom () {};
    Point_Custom (const double x, const double y, const double z)
    {
        this->x = x;
        this->y = y;
        this->z = z;
    };
    void print ()
    {
        printf("(%.2lf,%.2lf,%.2lf)\n",this->x,this->y,this->z);
    };
    bool operator <(const Point_Custom& pt) const
    {
        return (this->x < pt.x) || \
               ((!(pt.x < this->x)) && (this->y < pt.y)) || \
               ((!(pt.x < this->x)) && (!(pt.y < this->y)) && (this->z < pt.z));
    }
};

class Face_Custom
{
public:
    Point_Custom *v1;
    Point_Custom *v2;
    Point_Custom *v3;
    double normal[3];
public:
    Face_Custom (Point_Custom *a, Point_Custom *b, Point_Custom *c, double n[])
    {
        this->v1 = a;
        this->v2 = b;
        this->v3 = c;
        for (uint32_t i = 0; i < 3; i++)
            this->normal[i] = n[i];
    };
    void print ()
    {
        printf("\tVertex = [%lf %lf %lf] [%lf %lf %lf] [%lf %lf %lf]\n",this->v1->x,this->v1->y,this->v1->z,\
                                                        this->v2->x,this->v2->y,this->v2->z,\
                                                        this->v3->x,this->v3->y,this->v3->z);
        printf("\tNormal = (%lf %lf %lf)\n",this->normal[0],this->normal[1],this->normal[2]);
    };
};

// GLOBAL VARIABLES
bool first_call = true;
std::vector<Face_Custom> mesh_faces;
std::map<Point_Custom,uint32_t> unique_points;
std::vector< std::vector< uint32_t > > points_to_faces;

// Functions headers
extern "C" void move (struct walker *the_walker, struct user_options *the_options);
extern "C" void respawn (struct walker_config *the_walker_config, double pos[]);
extern "C" void draw (struct walker_config *the_walker_config, const double root_pos[]);

// Auxiliary functions
void read_faces_from_stl(const char filename[]);
void read_points_from_faces_to_map (const char filename[]);
void insert_points_from_faces_to_map ();
void read_face (FILE *file, std::vector<Face_Custom> &faces);
void print_faces (std::vector<Face_Custom> faces);

#endif