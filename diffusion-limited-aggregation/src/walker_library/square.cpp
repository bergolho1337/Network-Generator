#include "square.h"

SET_WALKER_MOVE_FUNCTION (move)
{
    //printf("[square] Inside square walker move function ...\n");

    double width;
    get_parameter_value_from_map(the_options->walker_config->params,"width",&width);

    double height;
    get_parameter_value_from_map(the_options->walker_config->params,"height",&height);

    double dx = generate_random_number();
    double dy = generate_random_number();

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

    uint8_t type = rand() % 4;
    double number = (double) rand() / (double) RAND_MAX;
    switch (type)
    {
        case 0: {
                    pos[0] = number * WIDTH;
                    pos[1] = 0.0;
                    pos[2] = 0.0;
                    break;
                }
        case 1: {
                    pos[0] = number * WIDTH;
                    pos[1] = HEIGHT;
                    pos[2] = 0.0;
                    break;
                }
        case 2: {
                    pos[0] = 0.0;
                    pos[1] = number * HEIGHT;
                    pos[2] = 0.0;
                    break;
                }
        case 3: {
                    pos[0] = WIDTH;
                    pos[1] = number * HEIGHT;
                    pos[2] = 0.0;
                    break;
                }
    }
}

bool is_inside (const double width, const double height, const double x, const double y)
{
    if ( (x >= 0.0 && x <= width) && (y >= 0.0 && y <= height) )
        return true;
    else
        return false;
}