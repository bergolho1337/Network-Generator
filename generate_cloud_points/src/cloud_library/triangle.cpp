#include "triangle.h"

SET_CLOUD_GENERATOR (default_cloud_triangle)
{
    printf("\n[triangle] Generating default triangle cloud of points\n");

    double area = generator->area;
    uint32_t num_points = generator->num_points;
    double side_length = sqrt( 4.0*area/sqrt(3.0) );
    std::vector<Point> *points = generator->points;

    printf("[triangle] Side length = %g\n\n",side_length);
    draw_triangle_area(side_length);

    double pos[3];
    for (uint32_t i = 0; i < num_points; i++)
    {
        do
        {
            generate_point_inside_triangle(pos,side_length);
        }while (!check_point(pos,*points,0.001));
        
        Point p(pos[0],pos[1],pos[2]);
        points->push_back(p);
        
        printf("[triangle] Generating point = %u\n",i);
    }

}

SET_CLOUD_GENERATOR (spaced_cloud_triangle)
{

    printf("\n[triangle] Generating sparse triangle cloud of points\n");

    double area = generator->area;
    uint32_t num_points = generator->num_points;
    double side_length = sqrt( 4.0*area/sqrt(3.0) );
    std::vector<Point> *points = generator->points;

    double tolerance;
    get_parameter_value_from_map(config->param,"tolerance",&tolerance);
    printf("[triangle] Tolerance = %g\n",tolerance);

    printf("[triangle] Side length = %g\n\n",side_length);
    draw_triangle_area(side_length);

    double pos[3];
    for (uint32_t i = 0; i < num_points; i++)
    {
        do
        {
            generate_point_inside_triangle(pos,side_length);
        }while (!check_point(pos,*points,tolerance));
        
        Point p(pos[0],pos[1],pos[2]);
        points->push_back(p);
        
        printf("[triangle] Generating point = %u\n",i);
    }
    
}

void generate_point_inside_triangle (double pos[], const double side_length)
{
    bool ok = false;
    while (!ok)
    {
        double teta = generate_random_number()*(M_PI/3.0) + (4.0*M_PI/3.0);
        double r = generate_random_number()*side_length;

        pos[0] = 0 + r*cos(teta);
        pos[1] = 0 + r*sin(teta);
        pos[2] = 0 + 0;

        if (pos[1] > -side_length*sin(M_PI/3.0))
            ok = true;
    }
}

void draw_triangle_area (const double side_length)
{
    double l = side_length;
    double l2 = side_length / 2.0;
    
    vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
    
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    points->InsertNextPoint(0,0,0);
    points->InsertNextPoint(l2,-l*sin(M_PI/3.0),0);
    points->InsertNextPoint(-l2,-l*sin(M_PI/3.0),0);
    polydata->SetPoints(points);

    vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();
    vtkSmartPointer<vtkLine> line1 = vtkSmartPointer<vtkLine>::New();
    line1->GetPointIds()->SetId(0,0);
    line1->GetPointIds()->SetId(1,1);
    vtkSmartPointer<vtkLine> line2 = vtkSmartPointer<vtkLine>::New();
    line2->GetPointIds()->SetId(0,1);
    line2->GetPointIds()->SetId(1,2);
    vtkSmartPointer<vtkLine> line3 = vtkSmartPointer<vtkLine>::New();
    line3->GetPointIds()->SetId(0,2);
    line3->GetPointIds()->SetId(1,0);
    lines->InsertNextCell(line1);
    lines->InsertNextCell(line2);
    lines->InsertNextCell(line3);
    polydata->SetLines(lines);

    vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
    writer->SetFileName("output/triangle_area_cloud.vtp");
    writer->SetInputData(polydata);
    writer->Write();

}