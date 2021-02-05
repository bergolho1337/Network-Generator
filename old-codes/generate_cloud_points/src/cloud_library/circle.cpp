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

SET_CLOUD_GENERATOR (default_circumscribed_cloud_circle)
{
    printf("\n[circle] Generating circumscribed circle cloud of points\n");

    double area = generator->area;
    uint32_t num_points = generator->num_points;
    double inner_radius = INNER_RADIUS;
    double inner_area = M_PI * inner_radius * inner_radius;
    double outer_radius = sqrt( (area + inner_area) / M_PI );
    std::vector<Point> *points = generator->points;

    printf("[circle] Inner radius = %g\n\n",inner_radius);
    printf("[circle] Outer radius = %g\n\n",outer_radius);
    draw_circumscribed_circle_area(inner_radius,outer_radius);

    double pos[3];
    for (uint32_t i = 0; i < num_points; i++)
    {
        do
        {
            generate_point_inside_circumscribed_circle(pos,inner_radius,outer_radius);
        }while (!check_point(pos,*points,0.001));
        
        Point p(pos[0],pos[1],pos[2]);
        points->push_back(p);
        
        printf("[circle] Generating point = %u\n",i);
    }
}

SET_CLOUD_GENERATOR (spaced_circumscribed_cloud_circle)
{
    printf("\n[circle] Generating spaced circumscribed circle cloud of points\n");

    double area = generator->area;
    uint32_t num_points = generator->num_points;
    double inner_radius = INNER_RADIUS;
    double inner_area = M_PI * inner_radius * inner_radius;
    double outer_radius = sqrt( (area + inner_area) / M_PI );
    std::vector<Point> *points = generator->points;

    double tolerance;
    get_parameter_value_from_map(config->param,"tolerance",&tolerance);
    printf("[circle] Tolerance = %g\n",tolerance);

    printf("[circle] Inner radius = %g\n\n",inner_radius);
    printf("[circle] Outer radius = %g\n\n",outer_radius);
    draw_circumscribed_circle_area(inner_radius,outer_radius);

    double pos[3];
    for (uint32_t i = 0; i < num_points; i++)
    {
        do
        {
            generate_point_inside_circumscribed_circle(pos,inner_radius,outer_radius);
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

    static double center[3] = {1,1,0};

    pos[0] = center[0] + r*cos(teta);
    pos[1] = center[1] + r*sin(teta);
    pos[2] = center[2] + 0;
}

void generate_point_inside_circumscribed_circle (double pos[],\
                                            const double inner_radius, const double outer_radius)
{

    bool ok = false;
    static const double center[3] = {0,-outer_radius,0};

    while (!ok)
    {
        double teta = generate_random_number()*2.0*M_PI;
        double r = generate_random_number()*outer_radius;

        pos[0] = 0 + r*cos(teta);
        pos[1] = -outer_radius + r*sin(teta);
        pos[2] = 0 + 0;

        double d = calc_euclidean_dist(pos,center);
        if (d >= inner_radius)
            ok = true;    
    }
    
}

void draw_circle_area (const double radius)
{
    vtkSmartPointer<vtkRegularPolygonSource> polygonSource =
      vtkSmartPointer<vtkRegularPolygonSource>::New();
    
    polygonSource->GeneratePolygonOff(); 
    polygonSource->SetNumberOfSides(50);
    polygonSource->SetRadius(radius);
    polygonSource->SetCenter(radius,radius,0);

    vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
    writer->SetFileName("output/circle_area_cloud.vtp");
    writer->SetInputConnection(polygonSource->GetOutputPort());
    writer->Write();
}

void draw_circumscribed_circle_area (const double inner_radius, const double outer_radius)
{
    vtkSmartPointer<vtkAppendPolyData> appendFilter = vtkSmartPointer<vtkAppendPolyData>::New();

    vtkSmartPointer<vtkRegularPolygonSource> inner_circle =
      vtkSmartPointer<vtkRegularPolygonSource>::New();
    
    inner_circle->GeneratePolygonOff(); 
    inner_circle->SetNumberOfSides(50);
    inner_circle->SetRadius(inner_radius);
    inner_circle->SetCenter(0,-outer_radius,0);

    vtkSmartPointer<vtkRegularPolygonSource> outer_circle =
      vtkSmartPointer<vtkRegularPolygonSource>::New();
    
    outer_circle->GeneratePolygonOff(); 
    outer_circle->SetNumberOfSides(50);
    outer_circle->SetRadius(outer_radius);
    outer_circle->SetCenter(0,-outer_radius,0);

    appendFilter->AddInputConnection(inner_circle->GetOutputPort());
    appendFilter->AddInputConnection(outer_circle->GetOutputPort());
    appendFilter->Update();

    vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
    writer->SetFileName("output/circumscribed_circle_area_cloud.vtp");
    writer->SetInputConnection(appendFilter->GetOutputPort());
    writer->Write();
}