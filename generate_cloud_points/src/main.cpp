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
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkDoubleArray.h>
#include <vtkAppendPolyData.h>
#include <vtkSphereSource.h>
#include <vtkCylinderSource.h>
#include <vtkMath.h>
#include <vtkLine.h>

using namespace std;

const unsigned int NUM_POINTS = 400;        // Total number of points to be generated
const double TOLERANCE = 0.3;               // Distance tolerance of the points
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

void draw_square_perfusion_area (const double length)
{

    // Insert the points
    double l = length;
    double l2 = length / 2.0;

    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    points->InsertNextPoint(-l2,l2,0);
    points->InsertNextPoint(l2,l2,0);
    points->InsertNextPoint(-l2,-l2,0);
    points->InsertNextPoint(l2,-l2,0);
    
    vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
    polydata->SetPoints(points);

    // Insert the lines using the edges of the graph
    vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();
    vtkSmartPointer<vtkLine> line1 = vtkSmartPointer<vtkLine>::New();
    line1->GetPointIds()->SetId(0,0);
    line1->GetPointIds()->SetId(1,1);

    vtkSmartPointer<vtkLine> line2 = vtkSmartPointer<vtkLine>::New();
    line2->GetPointIds()->SetId(0,1);
    line2->GetPointIds()->SetId(1,3);

    vtkSmartPointer<vtkLine> line3 = vtkSmartPointer<vtkLine>::New();
    line3->GetPointIds()->SetId(0,3);
    line3->GetPointIds()->SetId(1,2);

    vtkSmartPointer<vtkLine> line4 = vtkSmartPointer<vtkLine>::New();
    line4->GetPointIds()->SetId(0,2);
    line4->GetPointIds()->SetId(1,0);

    lines->InsertNextCell(line1);
    lines->InsertNextCell(line2);
    lines->InsertNextCell(line3);
    lines->InsertNextCell(line4);

    polydata->SetLines(lines);

    // Write the polydata to a file
    vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
    writer->SetFileName("output/perfusion_area_cloud.vtp");
    writer->SetInputData(polydata);
    writer->Write();
}

double generate_random_number ()
{
    // Between 0 and 1
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

void generate_point_inside_square (double pos[], const double length)
{
    int sign;
    double size = length / 2.0;

    double x = generate_random_number()*size;
    sign = rand() % 2;
    if (sign) 
        x *= -1.0;

    double y = generate_random_number()*size;
    sign = rand() % 2;
    if (sign) 
        y *= -1.0;

    pos[0] = x;
    pos[1] = y;
    pos[2] = 0.0;
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

void generate_square_cloud_points (vector<Point> &points, const double length)
{
    double pos[3];

    for (unsigned int i = 0; i < NUM_POINTS; i++)
    {
        do
        {
            generate_point_inside_square(pos,length);
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

void write_to_vtp (const vector<Point> points, const double radius)
{
    vtkSmartPointer<vtkAppendPolyData> appendFilter = vtkSmartPointer<vtkAppendPolyData>::New();
    
    // For each point create a sphere
    for (unsigned int i = 0; i < points.size(); i++)
    {
        vtkSmartPointer<vtkSphereSource> sphereSource = vtkSmartPointer<vtkSphereSource>::New();
        sphereSource->SetCenter(points[i].x,points[i].y,points[i].z);
        sphereSource->SetRadius(0.01*radius);
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

void generate_circle_cloud (const double A_perf)
{
    printf("[generate_cloud] Generating circle cloud of points\n");

    double r_perf = sqrt(A_perf / M_PI);

    printf("Radius = %g\n",r_perf);
    draw_perfusion_area(r_perf);

    vector<Point> points;
    generate_cloud_points(points,r_perf);

    write_to_vtp(points,r_perf);
    write_to_txt(points);
}

void generate_square_cloud (const double A_perf)
{
    printf("[generate_cloud] Generating square cloud of points\n");
    
    double l = sqrt(A_perf);
    printf("Square side = %g\n",l);

    draw_square_perfusion_area(l);

    vector<Point> points;
    generate_square_cloud_points(points,l);    

    write_to_vtp(points,l);
    write_to_txt(points);
}

void generate_triangle_cloud (const double A_perf)
{
    printf("[generate_cloud] Generating triangle cloud of points\n");
    
    printf("Need to implement !\n");
}

void generate_double_circle_cloud (const double A_perf)
{
    printf("[generate_cloud] Generating double circle cloud of points\n");
    
    printf("Need to implement !\n");
}

int main (int argc, char *argv[])
{
    if (argc-1 != 3)
    {
        printf("Usage:> %s <A_perf> <N_term> <type_cloud>\n",argv[0]);
        exit(EXIT_FAILURE);   
    }
    double A_perf = atof(argv[1]);
    int N_term = atoi(argv[2]);
    int type_cloud = atoi(argv[3]);

    switch (type_cloud)
    {
        case 0:
            generate_circle_cloud(A_perf);
            break;
        case 1:
            generate_square_cloud(A_perf);
            break;
        case 2:
            generate_triangle_cloud(A_perf);
            break;
        case 4:
            generate_double_circle_cloud(A_perf);
            break;

        default:
            fprintf(stderr,"[generate_cloud] Error! Invalid option !\n");
            exit(EXIT_FAILURE);
    }

    return 0;
}