#include "box.h"

SET_WALKER_MOVE_FUNCTION (move)
{
    double width = 600.0;
    get_parameter_value_from_map(the_options->walker_config->params,"width",&width);

    double height = 600.0;
    get_parameter_value_from_map(the_options->walker_config->params,"height",&height);

    double depth = 600.0;
    get_parameter_value_from_map(the_options->walker_config->params,"depth",&depth);

    double move_scale = 1.0;
    get_parameter_value_from_map(the_options->walker_config->params,"move_scale",&move_scale);

    double dx = generate_random_number()*move_scale;
    double dy = generate_random_number()*move_scale;
    double dz = generate_random_number()*move_scale;

    double new_x = the_walker->pos[0] + dx;
    double new_y = the_walker->pos[1] + dy;
    double new_z = the_walker->pos[2] + dz;

    if (is_inside(width,height,depth,new_x,new_y,new_z))
    {
        the_walker->pos[0] += dx;
        the_walker->pos[1] += dy;
        the_walker->pos[2] += dz;
    }
}

SET_WALKER_RESPAWN_FUNCTION (respawn)
{

    double width;
    get_parameter_value_from_map(the_walker_config->params,"width",&width);

    double height;
    get_parameter_value_from_map(the_walker_config->params,"height",&height);

    double depth = 600.0;
    get_parameter_value_from_map(the_walker_config->params,"depth",&depth);

    uint8_t type = rand() % 6;
    double number_1, number_2;
    number_1 = (double) rand() / (double) RAND_MAX;
    number_2 = (double) rand() / (double) RAND_MAX;
    switch (type)
    {
        case 0: {
                    pos[0] = number_1 * width;
                    pos[1] = 0.0;
                    pos[2] = number_2 * depth;
                    break;
                }
        case 1: {
                    pos[0] = 0.0;
                    pos[1] = number_1 * height;
                    pos[2] = number_2 * depth;
                    break;
                }
        case 2: {
                    pos[0] = number_1 * width;
                    pos[1] = number_2 * height;
                    pos[2] = depth;
                    break;
                }
        case 3: {
                    pos[0] = width;
                    pos[1] = number_1 * height;
                    pos[2] = number_2 * depth;
                    break;
                }
        case 4: {
                    pos[0] = number_1 * width;
                    pos[1] = number_2 * height;
                    pos[2] = 0.0;
                    break;
                }
        case 5: {
                    pos[0] = number_1 * width;
                    pos[1] = height;
                    pos[2] = number_2 * depth;
                    break;
                }
    }
}

SET_WALKER_DOMAIN_DRAW_FUNCTION (draw)
{
    double width = 600.0;
    get_parameter_value_from_map(the_walker_config->params,"width",&width);

    double height = 600.0;
    get_parameter_value_from_map(the_walker_config->params,"height",&height);

    double depth = 600.0;
    get_parameter_value_from_map(the_walker_config->params,"depth",&depth);

    double P0[3] = {0, 0, 0};
    double P1[3] = {0, 0, depth};
    double P2[3] = {width, 0, depth};
    double P3[3] = {width, 0, 0};
    double P4[3] = {0, height, 0};
    double P5[3] = {0, height, depth};
    double P6[3] = {width, height, depth};
    double P7[3] = {width, height, 0};

    // Create the points
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    points->InsertNextPoint(P0);
    points->InsertNextPoint(P1);
    points->InsertNextPoint(P2);
    points->InsertNextPoint(P3);
    points->InsertNextPoint(P4);
    points->InsertNextPoint(P5);
    points->InsertNextPoint(P6);
    points->InsertNextPoint(P7);

    // Create a hexahedron from the points
    vtkSmartPointer<vtkHexahedron> hex = vtkSmartPointer<vtkHexahedron>::New();
    hex->GetPointIds()->SetId(0,0);
    hex->GetPointIds()->SetId(1,1);
    hex->GetPointIds()->SetId(2,2);
    hex->GetPointIds()->SetId(3,3);
    hex->GetPointIds()->SetId(4,4);
    hex->GetPointIds()->SetId(5,5);
    hex->GetPointIds()->SetId(6,6);
    hex->GetPointIds()->SetId(7,7);

    // Add the hexahedron to a cell array
    vtkSmartPointer<vtkCellArray> hexs = vtkSmartPointer<vtkCellArray>::New();
    hexs->InsertNextCell(hex);

    // Add the points and hexahedron to an unstructured grid
    vtkSmartPointer<vtkUnstructuredGrid> unstructured_grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
    unstructured_grid->SetPoints(points);
    unstructured_grid->InsertNextCell(hex->GetCellType(), hex->GetPointIds());

    vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
    writer->SetFileName("output/box_walker_region.vtu");
    #if VTK_MAJOR_VERSION <= 5
    writer->SetInput(unstructured_grid);
    #else
    writer->SetInputData(unstructured_grid);
    #endif
    writer->Write();
}

bool is_inside (const double width, const double height, const double depth,\
                const double x, const double y, const double z)
{
    if ( (x >= 0.0 && x <= width) && (y >= 0.0 && y <= height) && (z >= 0.0 && z <= depth) )
        return true;
    else
        return false;
}