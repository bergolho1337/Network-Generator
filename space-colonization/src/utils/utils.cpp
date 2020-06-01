#include "utils.h"

void usage (const char pname[])
{
    fprintf(stderr,"%s\n",PRINT_LINE);
    fprintf(stderr,"Usage:> %s <input_config_file>\n",pname);
    fprintf(stderr,"%s\n",PRINT_DOTS);
    fprintf(stderr,"Examples:\n");
    fprintf(stderr,"\t%s inputs/simple_square_domain.ini\n",pname);
    fprintf(stderr,"\t%s inputs/simple_circle_domain.ini\n",pname);
    fprintf(stderr,"\t%s inputs/simple_box_domain.ini\n",pname);
    fprintf(stderr,"\t%s inputs/simple_sphere_domain.ini\n",pname);
    fprintf(stderr,"\t%s inputs/simple_cylinder_domain.ini\n",pname);
    fprintf(stderr,"\t%s inputs/simple_custom_domain.ini\n",pname);
    fprintf(stderr,"%s\n",PRINT_LINE);
}

void calculate_new_branch_position (double new_pos[], const double cur_pos[], const double u[], const double length, const double angle_degree)
{
    double angle_radians = convert_degrees_to_radians(angle_degree);

    // Rotation
    double r[3];
    r[0] = u[0]*cos(angle_radians) - u[1]*sin(angle_radians);
    r[1] = u[0]*sin(angle_radians) + u[1]*cos(angle_radians);
    r[2] = 0.0;

    new_pos[0] = cur_pos[0] + r[0]*length;
    new_pos[1] = cur_pos[1] + r[1]*length;
    new_pos[2] = cur_pos[2] + r[2]*length;
}