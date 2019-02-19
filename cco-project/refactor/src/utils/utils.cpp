#include "utils.h"

double euclidean_norm (const double x1, const double y1, const double z1,\
                    const double x2, const double y2, const double z2)
{
    return sqrt( pow(x2-x1,2) + pow(y2-y1,2) + pow(z2-z1,2) );
}

double generate_random_number ()
{
    double number = (double)rand() / (double)RAND_MAX;
    double sign = rand() % 2;
    
    return (sign) ? number : -number;
}

void generate_point_inside_perfusion_area (double pos[], const double radius)
{
    double teta = generate_random_number()*2.0*M_PI;
    double r = generate_random_number();

    // Center of perfusion area = (0,-radius,0)
    pos[0] = r*cos(teta);
    pos[1] = -r + r*sin(teta);
    pos[2] = 0;
}