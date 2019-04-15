// Author: Lucas Berg

#include <iostream>
#include <vector>
#include <random>
#include <algorithm>

#include <cstdio>
#include <cstdlib>
#include <cmath>

#include <vtkPoints.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkAppendPolyData.h>
#include <vtkSphereSource.h>

using namespace std;

// ====================================================================================================================
// CONSTANTS AND MACROS
#define PRINT_LINE "--------------------------------------------------------------------------------------------"

const double NE = 16;

class Point
{
public:
    double x, y, z;
public:
    Point (const double x, const double y, const double z)
    {
        this->x = x;
        this->y = y;
        this->z = z;
    }
    void print ()
    {
        printf("(%g,%g,%g)\n",this->x,this->y,this->z);
    }
};

// ====================================================================================================================

void create_sphere (vtkSmartPointer<vtkAppendPolyData> appendFilter, const double x[])
{
    vtkSmartPointer<vtkSphereSource> sphereSource = vtkSmartPointer<vtkSphereSource>::New();
    sphereSource->SetCenter(x[0],x[1],x[2]);
    sphereSource->SetRadius(0.05);
    sphereSource->Update();

    appendFilter->AddInputConnection(sphereSource->GetOutputPort());
    appendFilter->Update();
}

void read_user_input (vector<Point> &points, const char filename[])
{
    FILE *file = fopen(filename,"r");

    for (uint32_t i = 0; i < 3; i++)
    {
        double pos[3];
        fscanf(file,"%lf %lf %lf",&pos[0],&pos[1],&pos[2]);
        Point p(pos[0],pos[1],pos[2]);
        points.push_back(p);
    }

    fclose(file);
}

void print_points (vector<Point> points)
{
    printf("Number of points = %u\n",points.size());
    for (uint32_t i = 0; i < points.size(); i++)
        points[i].print();
}

bool is_corner (const uint32_t i, const uint32_t j)
{
    if ((i == 0 && j == 0) || (i == 0 && j == NE) || (i == NE && j == 0))
        return true;
    else
        return false;
}

void build_phi (double phi[], const double epsilon, const double eta)
{
    phi[0] = 1.0 - epsilon - eta;
    phi[1] = epsilon;
    phi[2] = eta;
}

void isoparametric_mapping (vector<Point> &mapped_points, vector<Point> main_points, const double delta_e)
{
    for (uint32_t i = 0; i <= NE; i++)
    {
        for (uint32_t j = 0; j <= NE-i; j++)
        {
            double epsilon = i*delta_e;
            double eta = j*delta_e;

            double phi[3];
            build_phi(phi,epsilon,eta);

            if (!is_corner(i,j))
            {
                double pos[3];
                pos[0] = (phi[0]*main_points[0].x) + (phi[1]*main_points[1].x) + (phi[2]*main_points[2].x);
                pos[1] = (phi[0]*main_points[0].y) + (phi[1]*main_points[1].y) + (phi[2]*main_points[2].y);
                pos[2] = (phi[0]*main_points[0].z) + (phi[1]*main_points[1].z) + (phi[2]*main_points[2].z);
            
                Point p(pos[0],pos[1],pos[2]);
                mapped_points.push_back(p);
            }
        }
    }
}

void write_points_to_vtu (vector<Point> main_points, vector<Point> mapped_points)
{
    vtkSmartPointer<vtkAppendPolyData> appendFilter = vtkSmartPointer<vtkAppendPolyData>::New();

    for (uint32_t i = 0; i < main_points.size(); i++)
    {
        double pos[3];
        pos[0] = main_points[i].x;
        pos[1] = main_points[i].y;
        pos[2] = main_points[i].z;

        create_sphere(appendFilter,pos);
    }

    for (uint32_t i = 0; i < mapped_points.size(); i++)
    {
        double pos[3];
        pos[0] = mapped_points[i].x;
        pos[1] = mapped_points[i].y;
        pos[2] = mapped_points[i].z;

        create_sphere(appendFilter,pos);
    }

    vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
    writer->SetFileName("output/points.vtp");
    writer->SetInputConnection(appendFilter->GetOutputPort());
    writer->Write();

}

int main (int argc, char *argv[])
{
    if (argc-1 != 1)
    {
        printf("%s\n",PRINT_LINE);
        printf("Usage:> %s <input_file>\n",argv[0]);
        printf("%s\n",PRINT_LINE);
        printf("<input_file> = Input file with the points coordinates\n");
        exit(EXIT_FAILURE);   
    }

    double delta_e = 1.0 / NE;

    vector<Point> main_points;
    read_user_input(main_points,argv[1]);
    //print_points(main_points);

    vector<Point> mapped_points;
    isoparametric_mapping(mapped_points,main_points,delta_e);
    //print_points(mapped_points);

    write_points_to_vtu(main_points,mapped_points);

    return 0;
}