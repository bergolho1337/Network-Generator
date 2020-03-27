#include "square.h"

// Generate points inside a square with dimensions: [-side_length,side_length]x[-side_length,side_length]x[0]
SET_CLOUD_GENERATOR (default_square_cloud)
{
    printf("\n[square] Generating default square cloud of points\n");

    uint32_t num_points = generator->num_points;
    std::vector<Point> *points = generator->points;
    double *rand_array = random_generator->array;
    
    double side_length;
    get_parameter_value_from_map(config->param,"side_length",&side_length);

    double pos[3];
    for (uint32_t i = 0; i < num_points; i++)
    {
        pos[0] = 2.0 * random_generator->get_value() - 1.0;
        pos[1] = 2.0 * random_generator->get_value() - 1.0;
        pos[2] = 0.0;

        // Convert to the real domain
        pos[0] *= side_length;
        pos[1] *= side_length;
        pos[2] *= side_length;

        Point p(pos[0],pos[1],pos[2]);
        points->push_back(p);
        
        printf("[square] Generating point = %u\n",i);
    }

    draw_default_square_volume(side_length);
}

void draw_default_square_volume (const double side_length)
{

}
