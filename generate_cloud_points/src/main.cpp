// Author: Lucas Berg
// Program that generates a cloud of points around a circle

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <vector>

#include <vtkRegularPolygonSource.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkSphereSource.h>
#include <vtkAppendPolyData.h>

using namespace std;

const unsigned int NUM_POINTS = 200;        // Total number of points to be generated
const double TOLERANCE = 0.05;               // Distance tolerance of the points
const unsigned int KTERM = 10;

class Point
{
public:
    double x, y, z;
public:
    Point () {};
    Point (const double x, const double y, const double z)
    {
        this->x = x;
        this->y = y;
        this->z = z;
    };
    void print ()
    {
        printf("(%.2lf,%.2lf,%.2lf)\n",this->x,this->y,this->z);
    };
};

void draw_perfusion_area (const double radius)
{
    vtkSmartPointer<vtkRegularPolygonSource> polygonSource =
      vtkSmartPointer<vtkRegularPolygonSource>::New();
    
    polygonSource->GeneratePolygonOff(); 
    polygonSource->SetNumberOfSides(50);
    polygonSource->SetRadius(radius);
    polygonSource->SetCenter(0,-radius,0);

    vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
    writer->SetFileName("output/perfusion_area_cloud.vtp");
    writer->SetInputConnection(polygonSource->GetOutputPort());
    writer->Write();
}

double generate_random_number ()
{
    double number = (double)rand() / (double)RAND_MAX;
    
    return number;
    
}

void generate_point_inside_circle (double pos[], const double radius)
{
    double teta = generate_random_number()*2.0*M_PI;
    double r = generate_random_number()*radius;

    pos[0] = 0 + r*cos(teta);
    pos[1] = -radius + r*sin(teta);
    pos[2] = 0 + 0;
}

double calc_euclidean_dist (const double a[], const double b[])
{
    return sqrt( pow(a[0]-b[0],2) + pow(a[1]-b[1],2) + pow(a[2]-b[2],2) );
}

bool check_point (const double new_pos[], vector<Point> points)
{
    for (int i = 0; i < (int)points.size(); i++)
    {
        double pos[3];
        pos[0] = points[i].x;
        pos[1] = points[i].y;
        pos[2] = points[i].z;

        double dist = calc_euclidean_dist(pos,new_pos);
        if (dist < TOLERANCE)
            return false;
    }
    return true;
}

void generate_cloud_points (vector<Point> &points, const double radius)
{
    double pos[3];

    for (unsigned int i = 0; i < NUM_POINTS; i++)
    {
        do
        {
            generate_point_inside_circle(pos,radius);
        }while (!check_point(pos,points));
        
        Point p(pos[0],pos[1],pos[2]);
        points.push_back(p);
        
        printf("[!] Number of points = %u\n",i);
    }    
}

void write_to_txt (const vector<Point> points)
{
    FILE *file = fopen("output/cloud_points.txt","w+");

    fprintf(file,"%lu\n",points.size());
    for (unsigned int i = 0; i < points.size(); i++)    
        fprintf(file,"%g %g %g\n",points[i].x,points[i].y,points[i].z);

    fclose(file);
}

void write_to_vtp (const vector<Point> points)
{
    vtkSmartPointer<vtkAppendPolyData> appendFilter = vtkSmartPointer<vtkAppendPolyData>::New();
    
    // For each point create a sphere
    for (unsigned int i = 0; i < points.size(); i++)
    {
        vtkSmartPointer<vtkSphereSource> sphereSource = vtkSmartPointer<vtkSphereSource>::New();
        sphereSource->SetCenter(points[i].x,points[i].y,points[i].z);
        sphereSource->SetRadius(0.01);
        sphereSource->Update();

        // Append the sphere to the filter
        appendFilter->AddInputConnection(sphereSource->GetOutputPort());
        appendFilter->Update();
    }

    // Write the file
    vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
    writer->SetFileName("output/cloud_points.vtp");
    writer->SetInputConnection(appendFilter->GetOutputPort());
    writer->Write();
}

int main (int argc, char *argv[])
{
    if (argc-1 != 2)
    {
        printf("Usage:> %s <A_perf> <N_term>\n",argv[0]);
        exit(EXIT_FAILURE);   
    }
    double A_perf = atof(argv[1]);
    int N_term = atoi(argv[2]);

    double A_sup = (double)((KTERM + 1) * A_perf) / (double)N_term;
    double r_sup = sqrt(A_sup / M_PI);

    printf("Radius = %g\n",r_sup);
    draw_perfusion_area(r_sup);

    vector<Point> points;
    generate_cloud_points(points,r_sup);

    write_to_vtp(points);
    write_to_txt(points);

    return 0;
}