#include "circle.h"

SET_WALKER_MOVE_FUNCTION (move)
{
    double *root_pos = the_options->root_pos;

    double radius = 300.0;
    get_parameter_value_from_map(the_options->walker_config->params,"radius",&radius);

    double move_scale = 5.0;
    get_parameter_value_from_map(the_options->walker_config->params,"move_scale",&move_scale);

    double dx = generate_random_number()*move_scale;
    double dy = generate_random_number()*move_scale;

    double new_pos[3];
    new_pos[0] = the_walker->pos[0] + dx;
    new_pos[1] = the_walker->pos[1] + dy;
    new_pos[2] = 0.0;

    if (is_inside(root_pos,radius,new_pos))
    {
        the_walker->pos[0] += dx;
        the_walker->pos[1] += dy;
    }
}

SET_WALKER_RESPAWN_FUNCTION (respawn)
{
    double radius = 300.0;
    get_parameter_value_from_map(the_walker_config->params,"radius",&radius);    

    double teta = ((double) rand() / (double) RAND_MAX)*2.0*M_PI;

    pos[0] = radius * cos(teta);
    pos[1] = radius * sin(teta);
    pos[2] = 0.0;
}

SET_WALKER_DOMAIN_DRAW_FUNCTION (draw)
{
    double radius;
    get_parameter_value_from_map(the_walker_config->params,"radius",&radius);
    
    vtkSmartPointer<vtkRegularPolygonSource> polygonSource = vtkSmartPointer<vtkRegularPolygonSource>::New();
    
    polygonSource->GeneratePolygonOff(); 
    polygonSource->SetNumberOfSides(50);
    polygonSource->SetRadius(radius);
    polygonSource->SetCenter(root_pos[0],root_pos[1],root_pos[2]);

    vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
    writer->SetFileName("outputs/circle_walker_region.vtp");
    writer->SetInputConnection(polygonSource->GetOutputPort());
    writer->Write();
}

bool is_inside (const double root_pos[], const double radius, const double new_pos[])
{
    double dist = calculate_distance(root_pos,new_pos);
    if (dist < radius * radius)
        return true;
    else
        return false;
}