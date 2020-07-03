// Author: Lucas Berg
// Program that receives a cloud of points and remap them using a root point as a reference.
// The remapping is done by envoloping the root point with a sphere, which at each iteration will
// increase its radius. By using this idea, all the current points of the cloud that are inside the spehre and 
// that are not visited yet will be placed in a array of points that can later be writen in .pts format.

#include <iostream>
#include <cmath>
#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <algorithm>

#include <vtkRegularPolygonSource.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkSphereSource.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkDoubleArray.h>
#include <vtkCylinderSource.h>
#include <vtkMath.h>
#include <vtkLine.h>
#include <vtkQuadric.h>
#include <vtkSampleFunction.h>
#include <vtkContourFilter.h>
#include <vtkPolyDataMapper.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkPolyDataWriter.h>
#include <vtkSTLWriter.h>
#include <vtkTransform.h>


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

void write_points_to_vtk (vector<Point> points, vector<bool> points_taken, const uint32_t iter)
{
    uint32_t num_points_taken = 0;
    for (uint32_t i = 0; i < points_taken.size(); i++)
    {
        if (!points_taken[i])
            num_points_taken++;
    }

    char filename[200];
    sprintf(filename,"outputs/points/points_%d.vtk",iter);
    FILE *file = fopen(filename,"w+");

    fprintf(file,"# vtk DataFile Version 3.0\n");
    fprintf(file,"Cloud\n");
    fprintf(file,"ASCII\n");
    fprintf(file,"DATASET POLYDATA\n");
    fprintf(file,"POINTS %lu float\n",num_points_taken);
    for (uint32_t i = 0; i < points_taken.size(); i++)
    {
        if (!points_taken[i])
            fprintf(file,"%lf %lf %lf\n",points[i].x,points[i].y,points[i].z);
    }
    fprintf(file,"VERTICES %lu %lu\n",num_points_taken,num_points_taken*2);
    for (uint32_t i = 0; i < num_points_taken; i++)
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

void draw_ellipsoid (const char filename[], const double center[], const double radius, const double a, const double b, const double c)
{
    // Create the quadric function definition
    vtkSmartPointer<vtkQuadric> quadric = vtkSmartPointer<vtkQuadric>::New();

    // Ellipsoid
    quadric->SetCoefficients(1/(a*a),1/(b*b),1/(c*c),0,0,0,0,0,0,-1);
    
    // Sample the quadric function
    vtkSmartPointer<vtkSampleFunction> sample = vtkSmartPointer<vtkSampleFunction>::New();
    sample->SetSampleDimensions(100,100,100);
    sample->SetImplicitFunction(quadric);
    double xmin = -5, xmax=5, ymin=-5, ymax=5, zmin=-5, zmax=5;
    //double xmin = -100, xmax=100, ymin=-100, ymax=100, zmin=-100, zmax=100;
    //xmin += new_radius; xmax += new_radius;
    //ymin += new_radius; ymax += new_radius;
    //zmin += new_radius; zmax += new_radius;
    //double xmin = -10, xmax=11, ymin=-10, ymax=10, zmin=-10, zmax=10;
    sample->SetModelBounds(xmin, xmax, ymin, ymax, zmin, zmax);

    //create the 0 isosurface (domain)
    vtkSmartPointer<vtkContourFilter> contours = vtkSmartPointer<vtkContourFilter>::New();
    contours->SetInputConnection(sample->GetOutputPort());
    contours->GenerateValues(1,1,1);

    vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->SetInputConnection(contours->GetOutputPort());
    mapper->Update();   

    // Get the reference to the Polydata
    vtkPolyData *polydata2 = mapper->GetInput();

    vtkSmartPointer<vtkTransform> translation = vtkSmartPointer<vtkTransform>::New();
    translation->Translate(center[0],center[1],center[2]);
    translation->Scale(radius*0.001,radius*0.001,radius*0.001);	

    vtkSmartPointer<vtkTransformPolyDataFilter> transformFilter = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
    transformFilter->SetInputData(polydata2);
    transformFilter->SetTransform(translation);
    transformFilter->Update();

    vtkSmartPointer<vtkPolyDataMapper> transformedMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    transformedMapper->SetInputConnection(transformFilter->GetOutputPort());
    transformedMapper->Update();

    // Get the reference to the Polydata
    //vtkPolyData *polydata = mapper->GetInput();
    vtkPolyData *polydata = transformedMapper->GetInput();

    // Write in VTK
    vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New();
    writer->SetFileName(filename);
    writer->SetInputData(polydata);
    writer->Write();

}

void calculate_points_limits (vector<Point> points, double &xmin, double &xmax,\
                             double &ymin, double &ymax, double &zmin, double &zmax)
{
    xmin = __DBL_MAX__; xmax = __DBL_MIN__;
    ymin = __DBL_MAX__; ymax = __DBL_MIN__;
    zmin = __DBL_MAX__; zmax = __DBL_MIN__;

    for (uint32_t i = 0; i < points.size(); i++)
    {
        if (points[i].x < xmin)
        {
            xmin = points[i].x;
        }
        if (points[i].x > xmax)
        {
            xmax = points[i].x;
        }

        if (points[i].y < ymin)
        {
            ymin = points[i].y;
        }
        if (points[i].y > ymax)
        {
            ymax = points[i].y;
        }

        if (points[i].z < zmin)
        {
            zmin = points[i].z;
        }
        if (points[i].z > zmax)
        {
            zmax = points[i].z;
        }
    }
}

void calculate_dimension_ratios (vector<Point> points, double ratios[])
{
    double xmin, xmax, ymin, ymax, zmin, zmax;
    double x_lim, y_lim, z_lim, min_lim, max_lim;
    double ratio_x, ratio_y, ratio_z;
    uint32_t max_axis, min_axis;
    
    calculate_points_limits(points,xmin, xmax, ymin, ymax, zmin, zmax);
    x_lim = fabs(xmin-xmax);
    y_lim = fabs(ymin-ymax);
    z_lim = fabs(zmin-zmax);
    max_lim = max(max(x_lim,y_lim),z_lim);

    ratio_x = max_lim/x_lim;
    ratio_y = max_lim/y_lim;
    ratio_z = max_lim/z_lim;
    
    ratios[0] = ratio_x;
    ratios[1] = ratio_y;
    ratios[2] = ratio_z;
    
    if (x_lim >= y_lim && x_lim >= z_lim)
        max_axis = 0;
    else if (y_lim >= x_lim && y_lim >= z_lim)
        max_axis = 1;
    else if (z_lim >= x_lim && z_lim >= y_lim)
        max_axis = 2;
    
    if (x_lim <= y_lim && x_lim <= z_lim)
        min_axis = 0;
    else if (y_lim <= x_lim && y_lim <= z_lim)
        min_axis = 1;
    else if (z_lim <= x_lim && z_lim <= y_lim)
        min_axis = 2;

    // Swap the minimum and maximum position
    double tmp;
    tmp = ratios[min_axis];
    ratios[min_axis] = ratios[max_axis];
    ratios[max_axis] = tmp;

    printf("------------------------------------------------\n");
    printf("xmin = %g || xmax = %g || x_lim = %g\n",xmin,xmax,x_lim);
    printf("ymin = %g || ymax = %g || y_lim = %g\n",ymin,ymax,y_lim);
    printf("zmin = %g || zmax = %g || z_lim = %g\n",zmin,zmax,z_lim);
    printf("------------------------------------------------\n");
    printf("ratio_x = %g\n",ratios[0]);
    printf("ratio_y = %g\n",ratios[1]);
    printf("ratio_z = %g\n",ratios[2]);
    printf("------------------------------------------------\n");
}

bool is_inside_ellipsoid (Point p, const double center[], const double radius, const double a, const double b, const double c)
{
    double r = radius*0.001;
    double dx = (p.x - center[0]);
    double dy = (p.y - center[1]);
    double dz = (p.z - center[2]);

    double value = ((dx*dx)/(a*a)) + ((dy*dy)/(b*b)) + ((dz*dz)/(c*c));
    if (value < r*r)
        return true;
    else
        return false;
}

void remap_points_inside_ellipsoid(const double center[], const double radius, const double a, const double b, const double c, vector<Point> &points, vector<Point> &remapped_points, vector<bool> &points_taken)
{
    uint32_t total_number_points = points.size(); 
    for (uint32_t i = 0; i < total_number_points; i++)
    {
        if (is_inside_ellipsoid(points[i],center,radius,a,b,c) && !points_taken[i])
        {
            remapped_points.push_back(points[i]);

            points_taken[i] = true;
            //points.erase(points.begin()+i);
        }
    }
}

void remap_points_from_root_v2 (vector<Point> &points, const uint32_t root_index)
{
    vector<Point> remapped_points;
    vector<bool> points_taken;

    double center[3];
    center[0] = points[root_index].x;
    center[1] = points[root_index].y;
    center[2] = points[root_index].z;

    points_taken.assign(points.size(),false);

    double ratios[3];
    calculate_dimension_ratios(points,ratios);

    double initial_radius = 1;
    double max_radius = 200;
    double offset = 2;
    //double a = ratios[0];         // X axis: || > 1 (decrease) || < 1 (increase)
    //double b = ratios[1];         // Y axis: || > 1 (decrease) || < 1 (increase)
    //double c = ratios[2];         // Z axis: || > 1 (decrease) || < 1 (increase)
    double a = 1.0;  // LV rabbit
    double b = 1.0;  // LV rabbit
    double c = 2.0;  // LV rabbit
    //double a = 0.5;  // RV rabbit
    //double b = 1.0;  // RV rabbit
    //double c = 5.0;  // RV rabbit
    uint32_t iter = 0;

    for (double radius = initial_radius; radius < max_radius; radius += offset, iter++)
    {
        printf("[remapping] Radius = %g ...\n",radius);

        char filename[200];
        sprintf(filename,"outputs/ellipsoids/ellipsoid_%d.vtk",iter);
        draw_ellipsoid(filename,center,radius,a,b,c);

        remap_points_inside_ellipsoid(center,radius,a,b,c,points,remapped_points,points_taken);
        
        write_points_to_vtk(points,points_taken,iter);
        write_remapped_points_to_vtk(remapped_points,iter);

    }
    write_points_to_pts(remapped_points);
}

int main (int argc, char *argv[])
{
    if (argc-1 != 2)
    {
        printf("----------------------------------------------------------------------------\n");
        printf("Usage:> %s <input_filename> <root_point_index>\n",argv[0]);
        printf("----------------------------------------------------------------------------\n");
        printf("Example:\n");
        printf("\t%s inputs/paraboloid_exterior_cloud_points.pts 90000\n",argv[0]);
        printf("\t%s inputs/elizabeth_exterior_LV.pts 1200\n",argv[0]);
        printf("----------------------------------------------------------------------------\n");
        exit(EXIT_FAILURE);   
    }

    vector<Point> points;
    read_points(argv[1],points);
    //print_points(points);

    uint32_t root_index = atoi(argv[2]);

    remap_points_from_root_v2(points,root_index);


    return 0;
}
