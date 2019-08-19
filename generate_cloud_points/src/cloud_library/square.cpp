#include "square.h"

SET_CLOUD_GENERATOR (default_cloud_square)
{
    printf("\n[square] Generating default square cloud of points\n");

    double area = generator->area;
    uint32_t num_points = generator->num_points;
    double side_length = sqrt(area);
    std::vector<Point> *points = generator->points;

    printf("[square] Side length = %g\n\n",side_length);
    draw_square_area(side_length);

    double pos[3];
    for (uint32_t i = 0; i < num_points; i++)
    {
        do
        {
            generate_point_inside_square(pos,side_length);
        }while (!check_point(pos,*points,0.001));
        
        Point p(pos[0],pos[1],pos[2]);
        points->push_back(p);
        
        printf("[square] Generating point = %u\n",i);
    }

}

SET_CLOUD_GENERATOR (spaced_cloud_square)
{

    printf("\n[square] Generating sparse square cloud of points\n");

    double area = generator->area;
    uint32_t num_points = generator->num_points;
    double side_length = sqrt(area);
    std::vector<Point> *points = generator->points;

    double tolerance;
    get_parameter_value_from_map(config->param,"tolerance",&tolerance);
    printf("[square] Tolerance = %g\n",tolerance);

    printf("Side length = %g\n\n",side_length);
    draw_square_area(side_length);

    double pos[3];
    for (uint32_t i = 0; i < num_points; i++)
    {
        do
        {
            generate_point_inside_square(pos,side_length);
        }while (!check_point(pos,*points,tolerance));
        
        Point p(pos[0],pos[1],pos[2]);
        points->push_back(p);
        
        printf("[square] Generating point = %u\n",i);
    }

}

void generate_point_inside_square (double pos[], const double side_length)
{
    double l2 = side_length / 2.0;

    double x = generate_random_number()*side_length;
    double y = generate_random_number()*side_length;

    pos[0] = x;
    pos[1] = y;
    pos[2] = 0.0;
}

void draw_square_area (const double side_length)
{
    double l = side_length;
    double l2 = side_length / 2.0;
    
    vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
    
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    points->InsertNextPoint(0,l,0);
    points->InsertNextPoint(l,l,0);
    points->InsertNextPoint(l,0,0);
    points->InsertNextPoint(0,0,0);
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
    line3->GetPointIds()->SetId(1,3);
    vtkSmartPointer<vtkLine> line4 = vtkSmartPointer<vtkLine>::New();
    line4->GetPointIds()->SetId(0,3);
    line4->GetPointIds()->SetId(1,0);
    lines->InsertNextCell(line1);
    lines->InsertNextCell(line2);
    lines->InsertNextCell(line3);
    lines->InsertNextCell(line4);
    polydata->SetLines(lines);

    vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
    writer->SetFileName("output/square_area_cloud.vtp");
    writer->SetInputData(polydata);
    writer->Write();

}
