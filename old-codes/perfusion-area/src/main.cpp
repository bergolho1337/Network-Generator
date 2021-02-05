// Author: Lucas Berg
// Program that draws the perfusion area

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <vector>

#include <vtkRegularPolygonSource.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkSphereSource.h>
#include <vtkAppendPolyData.h>

using namespace std;

// =============================================
// CONSTANTS AND MACROS
static const int NTERM = 10;
static const double APERF = 1.0;
// =============================================

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

void draw_perfusion_area (const char filename[], const double radius)
{
    vtkSmartPointer<vtkRegularPolygonSource> polygonSource =
      vtkSmartPointer<vtkRegularPolygonSource>::New();
    
    polygonSource->GeneratePolygonOff(); 
    polygonSource->SetNumberOfSides(50);
    polygonSource->SetRadius(radius);
    polygonSource->SetCenter(0,-radius,0);

    vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
    writer->SetFileName(filename);
    writer->SetInputConnection(polygonSource->GetOutputPort());
    writer->Write();
}

int main (int argc, char *argv[])
{
    if (argc-1 != 0)
    {
        printf("Usage:> %s\n",argv[0]);
        exit(EXIT_FAILURE);   
    }

    char filename[50];
    for (int Kterm = 1; Kterm < NTERM; Kterm++)
    {
        double radius = (double)((Kterm + 1) * APERF) / (double)NTERM;

        printf("Radius = %g\n",radius);

        sprintf(filename,"output/perfusion-area-%d.vtp",Kterm);

        draw_perfusion_area(filename,radius);
    }

    return 0;
}