#include "square.h"

SET_WALKER_MOVE_FUNCTION (move)
{
    //printf("[square] Inside square walker move function ...\n");

    double width = 600.0;
    get_parameter_value_from_map(the_options->walker_config->params,"width",&width);

    double height = 600.0;
    get_parameter_value_from_map(the_options->walker_config->params,"height",&height);

    double move_scale = 1.0;
    get_parameter_value_from_map(the_options->walker_config->params,"move_scale",&move_scale);

    double dx = generate_random_number()*move_scale;
    double dy = generate_random_number()*move_scale;

    double new_x = the_walker->pos[0] + dx;
    double new_y = the_walker->pos[1] + dy;

    if (is_inside(width,height,new_x,new_y))
    {
        the_walker->pos[0] += dx;
        the_walker->pos[1] += dy;
    }
}

SET_WALKER_RESPAWN_FUNCTION (respawn)
{
    //printf("[square] Inside square walker respawn function ...\n");

    double width;
    get_parameter_value_from_map(the_walker_config->params,"width",&width);

    double height;
    get_parameter_value_from_map(the_walker_config->params,"height",&height);

    uint8_t type = rand() % 4;
    double number = (double) rand() / (double) RAND_MAX;
    switch (type)
    {
        case 0: {
                    pos[0] = number * width;
                    pos[1] = 0.0;
                    pos[2] = 0.0;
                    break;
                }
        case 1: {
                    pos[0] = number * width;
                    pos[1] = height;
                    pos[2] = 0.0;
                    break;
                }
        case 2: {
                    pos[0] = 0.0;
                    pos[1] = number * height;
                    pos[2] = 0.0;
                    break;
                }
        case 3: {
                    pos[0] = width;
                    pos[1] = number * height;
                    pos[2] = 0.0;
                    break;
                }
    }
}

SET_WALKER_DOMAIN_DRAW_FUNCTION (draw)
{
    double width;
    get_parameter_value_from_map(the_walker_config->params,"width",&width);

    double height;
    get_parameter_value_from_map(the_walker_config->params,"height",&height);

    vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
    
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    points->InsertNextPoint(0,height,0);
    points->InsertNextPoint(width,height,0);
    points->InsertNextPoint(width,0,0);
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
    writer->SetFileName("output/square_walker_region.vtp");
    writer->SetInputData(polydata);
    writer->Write();
}

bool is_inside (const double width, const double height, const double x, const double y)
{
    if ( (x >= 0.0 && x <= width) && (y >= 0.0 && y <= height) )
        return true;
    else
        return false;
}