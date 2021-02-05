#include "utils.h"

double calc_norm (const double x1, const double y1, const double z1,\
                const double x2, const double y2, const double z2)
{
    return sqrt( pow(x1-x2,2) + pow(y1-y2,2) + pow(z1-z2,2) );
}

void build_unitary_vector (double d[], const double src[], const double dest[])
{
    double norm = calc_norm(src[0],src[1],src[2],\
                        dest[0],dest[1],dest[2]);
    if (norm < 1.0e-08)
        norm = 1.0e-08;

    for (uint32_t i = 0; i < 3; i++)
    {
        d[i] = (dest[i] - src[i]) / norm;
    }
        
}

double calc_angle_between_vectors (const double u[], const double v[])
{
    double dot_product = u[0]*v[0] + u[1]*v[1] + u[2]*v[2];

    double angle_radians = acos(dot_product);

    // Return the angle in degrees
    return angle_radians * 180.0 / M_PI;
}

void write_data_to_file (const char filename[], std::vector<double> arr)
{
    FILE *file = fopen(filename,"w+");

    for (uint32_t i = 0; i < arr.size(); i++)
        fprintf(file,"%g\n",arr[i]);

    fclose(file);
}

uint32_t calc_closest_point (Point p, std::vector<Point> points)
{
    double min_dist = __DBL_MAX__;
    uint32_t min_index = points.size()+1;
    for (uint32_t i = 0; i < points.size(); i++)
    {
        double dist = calc_norm(p.x,p.y,p.z,points[i].x,points[i].y,points[i].z);

        if (dist < min_dist)
        {
            min_dist = dist;
            min_index = i;
        }
    }
    return min_index;
}