#include "circle.h"

SET_CLOUD_GENERATOR (default_cloud_circle)
{
    printf("\n[circle] Generating default circle cloud of points\n");

    double area = generator->area;
    uint32_t num_points = generator->num_points;
    double radius = sqrt(area / M_PI);
    std::vector<Point> *points = generator->points;

    printf("Radius = %g\n\n",radius);
    draw_circle_area(radius);

    double pos[3];
    for (uint32_t i = 0; i < num_points; i++)
    {
        do
        {
            generate_point_inside_circle(pos,radius);
        }while (!check_point(pos,*points,0.001));
        
        Point p(pos[0],pos[1],pos[2]);
        points->push_back(p);
        
        printf("[circle] Generating point = %u\n",i);
    }

}

SET_CLOUD_GENERATOR (spaced_cloud_circle)
{

    printf("\n[circle] Generating sparse circle cloud of points\n");

    double area = generator->area;
    uint32_t num_points = generator->num_points;
    double radius = sqrt(area / M_PI);
    std::vector<Point> *points = generator->points;

    double tolerance;
    get_parameter_value_from_map(config->param,"tolerance",&tolerance);
    printf("[circle] Tolerance = %g\n",tolerance);

    printf("[circle] Radius = %g\n\n",radius);
    draw_circle_area(radius);

    double pos[3];
    for (uint32_t i = 0; i < num_points; i++)
    {
        do
        {
            generate_point_inside_circle(pos,radius);
        }while (!check_point(pos,*points,tolerance));
        
        Point p(pos[0],pos[1],pos[2]);
        points->push_back(p);
        
        printf("[circle] Generating point = %u\n",i);
    }
}

void generate_point_inside_circle (double pos[], const double radius)
{
    double teta = generate_random_number()*2.0*M_PI;
    double r = generate_random_number()*radius;

    pos[0] = 0 + r*cos(teta);
    pos[1] = -radius + r*sin(teta);
    pos[2] = 0 + 0;
}

void draw_circle_area (const double radius)
{
    vtkSmartPointer<vtkRegularPolygonSource> polygonSource =
      vtkSmartPointer<vtkRegularPolygonSource>::New();
    
    polygonSource->GeneratePolygonOff(); 
    polygonSource->SetNumberOfSides(50);
    polygonSource->SetRadius(radius);
    polygonSource->SetCenter(0,-radius,0);

    vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
    writer->SetFileName("output/circle_area_cloud.vtp");
    writer->SetInputConnection(polygonSource->GetOutputPort());
    writer->Write();
}