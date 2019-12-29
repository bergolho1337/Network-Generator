#include "utils.h"

void usage (const char pname[])
{
    fprintf(stderr,"%s\n",PRINT_LINE);
    fprintf(stderr,"Usage:> %s <input_config_file>\n",pname);
    fprintf(stderr,"%s\n",PRINT_LINE_2);
    fprintf(stderr,"<input_config_file> = Input file with user configurations\n");
    fprintf(stderr,"%s\n",PRINT_LINE);
}

double generate_random_number ()
{
    // Generate a number between 0 and 1
    return (double)rand() / (double)RAND_MAX;
}

bool check_point (const double new_pos[], std::vector<Point> points, const double tolerance)
{
    for (uint32_t i = 0; i < points.size(); i++)
    {
        double pos[3];
        pos[0] = points[i].x;
        pos[1] = points[i].y;
        pos[2] = points[i].z;

        double dist = calc_euclidean_dist(pos,new_pos);
        if (dist < tolerance)
            return false;
    }
    return true;
}

bool get_parameter_value_from_map (std::map<std::string,double> *params,\
                                const std::string key, double *value)
{
    auto it = params->find(key);

    if (it != params->end())
    {
        //printf("Found \"%s\" on parameter list: %g\n",key.c_str(),it->second);
        *value = it->second;
        return true;
    }
    else
    {
        fprintf(stderr,"Not found \"%s\" \n",key.c_str());
        return false;
    }
}

double calc_euclidean_dist (const double a[], const double b[])
{
    return sqrt( pow(a[0]-b[0],2) + pow(a[1]-b[1],2) + pow(a[2]-b[2],2) );
}