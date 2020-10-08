#include <iostream>
#include <cstdio>
#include <cmath>
#include <cstring>
#include <cstdlib>

#include "stl-reader/reader.h"

using namespace std;

void read_pmjs (const char filename[], vector<Point> &points)
{
    char str[200];
    int np;
    FILE *file = fopen(filename,"r");
    while (fscanf(file,"%s",str) != EOF)
    {
        if (strcmp(str,"POINTS") == 0) break;
    }
    fscanf(file,"%d %s",&np,str);
    for (int i = 0; i < np; i++)
    {
        double pos[3];
        fscanf(file,"%lf %lf %lf",&pos[0],&pos[1],&pos[2]);
        Point p(i,pos[0],pos[1],pos[2]);
        points.push_back(p);
    }
    fclose(file);
}

uint32_t find_closest_point (Point p, vector<Point> points)
{
    uint32_t min_id = 0;
    double min_dist = __DBL_MAX__;
    for (uint32_t i = 0; i < points.size(); i++)
    {
        double dist = sqrt(pow(p.x-points[i].x,2)+pow(p.y-points[i].y,2)+pow(p.z-points[i].z,2));
        if (dist < min_dist)
        {
            min_dist = dist;
            min_id = i;
        }
    }
    return min_id;
}

int main (int argc, char *argv[])
{
    if (argc-1 != 1)
    {
        printf("------------------------------------------------------------------------------------------------------------------------------------------\n");
        printf("Usage:> %s <STL_filename>\n",argv[0]);
        printf("------------------------------------------------------------------------------------------------------------------------------------------\n");
        exit(EXIT_FAILURE);   
    }

    // Read the STL mesh file and store its faces in a vector
    vector<Point> points;
    vector<Face> faces;
    read_faces(argv[1],points,faces);
    
    // This will print the STL file in the format that the Geodesic Library understand 
    printf("%u %u\n",points.size(),faces.size());
    print_points(points);
    print_faces(faces);

    vector<Point> pmj_points;
    read_pmjs("inputs/elizabeth_pmjs_RV.vtk",pmj_points);    
    
    printf("TARGETS=( ");
    for (uint32_t i = 0; i < pmj_points.size(); i++)
    {
        uint32_t id = find_closest_point(pmj_points[i],points);
        printf("%u ",id);
    }
    printf(")\n");

    return 0;
}
