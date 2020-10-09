#include <iostream>
#include <vector>
#include <algorithm>
#include <cstdio>
#include <cmath>
#include <cstring>
#include <cstdlib>

#include <vtkXMLPolyDataReader.h>
#include <vtkSTLReader.h>
#include <vtkPolyDataMapper.h>
#include <vtkSphereSource.h>
#include <vtkProperty.h>
#include <vtkPolyData.h>
#include <vtkIdList.h>
#include <vtkSmartPointer.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkDijkstraGraphGeodesicPath.h>

using namespace std;

const double SCALE_FACTOR = 0.25;

class Point
{
public:
    uint32_t id;
    double x, y, z;
public:
    Point (const uint32_t id, const double x, const double y, const double z)
    {
        this->id = id;
        this->x = x;
        this->y = y;
        this->z = z;
    }
};

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
    // Read the STL mesh file and store its points in a vector
    vtkSmartPointer<vtkSTLReader> reader = vtkSmartPointer<vtkSTLReader>::New();
    reader->SetFileName("inputs/elizabeth_endo_LV.stl");
    reader->Update();

    vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->SetInputConnection(reader->GetOutputPort());
    mapper->Update();

    vtkPolyData *polydata_grid = mapper->GetInput();

    vector<Point> points;
    uint32_t num_points = polydata_grid->GetNumberOfPoints();
    for (uint32_t i = 0; i < num_points; i++)
    {
        double pos[3];
        polydata_grid->GetPoint(i,pos);
        
        Point p(i,pos[0],pos[1],pos[2]);
        points.push_back(p);
    }

    // Read the PMJ locations file
    vector<Point> pmj_points;
    read_pmjs("inputs/elizabeth_pmjs_LV.vtk",pmj_points);    
    
    // Calculate the geodesic path from the root to each terminal
    vtkSmartPointer<vtkDijkstraGraphGeodesicPath> dijkstra = vtkSmartPointer<vtkDijkstraGraphGeodesicPath>::New();
    dijkstra->SetInputData(polydata_grid);
    dijkstra->SetStartVertex(753);      // LV
    //dijkstra->SetStartVertex(18325);    // RV
    vector<Point> path_points[pmj_points.size()];
    for (uint32_t i = 0; i < pmj_points.size(); i++)
    {
        uint32_t target = find_closest_point(pmj_points[i],points);
        
        dijkstra->SetEndVertex(target);
        dijkstra->Update();

        vtkSmartPointer<vtkIdList> ids = dijkstra->GetIdList();
        uint32_t num_ids = ids->GetNumberOfIds();
        
        for (uint32_t j = 0; j < num_ids; j++)
        {
            uint32_t id = ids->GetId(j);
            double pos[3];
            polydata_grid->GetPoint(id,pos);

            //printf("%u %g %g %g\n",id,pos[0],pos[1],pos[2]);
            Point p(j,pos[0],pos[1],pos[2]);
            path_points[i].push_back(p);
        }
        /*
        char filename[200];
        sprintf(filename,"pathways/terminal_%u.vtk",i);
        FILE *file = fopen(filename,"w+");
        fprintf(file,"# vtk DataFile Version 3.0\n");
        fprintf(file,"Geodesic Pathway\n");
        fprintf(file,"ASCII\n");
        fprintf(file,"DATASET POLYDATA\n");
        fprintf(file,"POINTS %u float\n",path_points[i].size());
        for (uint32_t j = 0; j < path_points[i].size(); j++)
            fprintf(file,"%g %g %g\n",path_points[i][j].x*SCALE_FACTOR,path_points[i][j].y*SCALE_FACTOR,path_points[i][j].z*SCALE_FACTOR);
        fprintf(file,"LINES %u %u\n",path_points[i].size()-1,(path_points[i].size()-1)*3);
        for (uint32_t j = 0; j < path_points[i].size()-1; j++)
            fprintf(file,"2 %u %u\n",j,j+1);
        fclose(file);
        */
    }

    char filename[200];
    sprintf(filename,"pathways/pathways.vtk");
    FILE *file = fopen(filename,"w+");
    
    uint32_t counter = 0;
    vector<Point> pathway_points;
    vector< pair<uint32_t,uint32_t> > pathway_lines;

    // Add the root
    uint32_t last = path_points[0].size()-1;
    Point p(counter,path_points[0][last].x,path_points[0][last].y,path_points[0][last].z);
    pathway_points.push_back(p);
    //pathway_lines.push_back(make_pair(counter,counter+1));
    counter++;

    for (uint32_t i = 0; i < pmj_points.size(); i++)
    {
        uint32_t first_index = counter;
        for (uint32_t j = path_points[i].size()-2; j > 0; j--)
        {
            Point p(counter,path_points[i][j].x,path_points[i][j].y,path_points[i][j].z);
            pathway_points.push_back(p);

            if (j > 1)
                pathway_lines.push_back(make_pair(counter,counter+1));

            counter++;
        }
        //printf("0 %u\n",first_index);
        pathway_lines.push_back(make_pair(0,first_index));
        
    }

    fprintf(file,"# vtk DataFile Version 3.0\n");
    fprintf(file,"Geodesic Pathway\n");
    fprintf(file,"ASCII\n");
    fprintf(file,"DATASET POLYDATA\n");
    fprintf(file,"POINTS %u float\n",pathway_points.size());
    for (uint32_t j = 0; j < pathway_points.size(); j++)
        fprintf(file,"%g %g %g\n",pathway_points[j].x*SCALE_FACTOR,pathway_points[j].y*SCALE_FACTOR,pathway_points[j].z*SCALE_FACTOR);
    fprintf(file,"LINES %u %u\n",pathway_lines.size(),pathway_lines.size()*3);
    for (uint32_t j = 0; j < pathway_lines.size(); j++)
        fprintf(file,"2 %u %u\n",pathway_lines[j].first,pathway_lines[j].second);
    fclose(file);
    
    return 0;
}
