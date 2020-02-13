// Author: Lucas Berg

#include <iostream>
#include <cmath>
#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <algorithm>

#include <vtkRegularPolygonSource.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkSphereSource.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkDoubleArray.h>
#include <vtkCylinderSource.h>
#include <vtkMath.h>
#include <vtkLine.h>


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

    fscanf(file,"%u",&num_points);
    
    for (uint32_t i = 0; i < num_points; i++)
    {
        fscanf(file,"%lf %lf %lf",&x,&y,&z);

        Point point(i,x,y,z);
        points.push_back(point);
    }

    fclose(file);
}

void print_points (vector<Point> points)
{
    for (int i = 0; i < (int)points.size(); i++)
        points[i].print();
}

void write_points_to_pts (vector<Point> points)
{
    char filename[200];
    sprintf(filename,"outputs/remapped_points.pts");
    FILE *file = fopen(filename,"w+");

    fprintf(file,"%lu\n",points.size());
    for (uint32_t i = 0; i < points.size(); i++)
        fprintf(file,"%lf %lf %lf\n",points[i].x,points[i].y,points[i].z);
    
    fclose(file);
}

void write_points_to_vtk (vector<Point> points, const uint32_t iter)
{
    char filename[200];
    sprintf(filename,"outputs/points/points_%d.vtk",iter);
    FILE *file = fopen(filename,"w+");

    fprintf(file,"# vtk DataFile Version 3.0\n");
    fprintf(file,"Cloud\n");
    fprintf(file,"ASCII\n");
    fprintf(file,"DATASET POLYDATA\n");
    fprintf(file,"POINTS %lu float\n",points.size());
    for (uint32_t i = 0; i < points.size(); i++)
        fprintf(file,"%lf %lf %lf\n",points[i].x,points[i].y,points[i].z);
    fprintf(file,"VERTICES %lu %lu\n",points.size(),points.size()*2);
    for (uint32_t i = 0; i < points.size(); i++)
        fprintf(file,"1 %u\n",i);
    
    fclose(file);
}

void write_remapped_points_to_vtk (vector<Point> points, const uint32_t iter)
{
    char filename[200];
    sprintf(filename,"outputs/remapped-points/remapped_points_%d.vtk",iter);
    FILE *file = fopen(filename,"w+");

    fprintf(file,"# vtk DataFile Version 3.0\n");
    fprintf(file,"Cloud\n");
    fprintf(file,"ASCII\n");
    fprintf(file,"DATASET POLYDATA\n");
    fprintf(file,"POINTS %lu float\n",points.size());
    for (uint32_t i = 0; i < points.size(); i++)
        fprintf(file,"%lf %lf %lf\n",points[i].x,points[i].y,points[i].z);
    fprintf(file,"VERTICES %lu %lu\n",points.size(),points.size()*2);
    for (uint32_t i = 0; i < points.size(); i++)
        fprintf(file,"1 %u\n",i);
    
    fclose(file);
}

void draw_sphere (const char filename[], const double center[], const double radius)
{
    vtkSmartPointer<vtkSphereSource> sphereSource = vtkSmartPointer<vtkSphereSource>::New();
    
    sphereSource->SetRadius(radius);
    sphereSource->SetCenter(center[0],center[1],center[2]);
    sphereSource->SetPhiResolution(20);
    sphereSource->SetThetaResolution(20);

    vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
    writer->SetFileName(filename);
    writer->SetInputConnection(sphereSource->GetOutputPort());
    writer->Write();
}

bool is_inside_sphere (Point p, const double center[], const double radius)
{
    double dx = (p.x - center[0]);
    double dy = (p.y - center[1]);
    double dz = (p.z - center[2]);

    double dist = sqrt((dx*dx) + (dy*dy) + (dz*dz));
    if (dist < radius)
        return true;
    else
        return false;
}

void remap_points_inside_sphere(const double center[], const double radius, vector<Point> &points, vector<Point> &remapped_points, vector<bool> &points_taken)
{
    uint32_t total_number_points = points.size(); 
    for (uint32_t i = 0; i < total_number_points; i++)
    {
        if (is_inside_sphere(points[i],center,radius) && !points_taken[i])
        {
            remapped_points.push_back(points[i]);

            points_taken[i] = true;
            //points.erase(points.begin()+i);
        }
    }
} 

void remap_points_from_root(vector<Point> &points, const uint32_t root_index)
{
    vector<Point> remapped_points;
    vector<bool> points_taken;

    double center[3];
    center[0] = points[root_index].x;
    center[1] = points[root_index].y;
    center[2] = points[root_index].z;

    points_taken.assign(points.size(),false);

    double initial_radius = 0.0001;
    double max_radius = 0.1;
    double offset = 0.001;
    uint32_t iter = 0;
    for (double radius = initial_radius; radius < max_radius; radius += offset, iter++)
    {
        char filename[200];
        sprintf(filename,"outputs/spheres/sphere_%d.vtp",iter);
        draw_sphere(filename,center,radius);

        remap_points_inside_sphere(center,radius,points,remapped_points,points_taken);
        
        write_points_to_vtk(points,iter);
        write_remapped_points_to_vtk(remapped_points,iter);

    }

    write_points_to_pts(remapped_points);
}

int main (int argc, char *argv[])
{
    if (argc-1 != 2)
    {
        printf("Usage:> %s <input_filename> <root_point_index>\n",argv[0]);
        exit(EXIT_FAILURE);   
    }

    vector<Point> points;
    read_points(argv[1],points);
    //print_points(points);

    uint32_t root_index = atoi(argv[2]);

    remap_points_from_root(points,root_index);


    return 0;
}
