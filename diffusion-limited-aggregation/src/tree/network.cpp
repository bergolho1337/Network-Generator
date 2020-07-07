#include "network.h"

void read_initial_network_file (const char filename[], std::vector<struct network_point> &points, std::vector<struct network_line> &lines)
{
    FILE *file = fopen(filename,"r");
    if (!file)
    {
        fprintf(stderr,"[-] ERROR! Reading initial network file '%s'\n",filename);
        exit(EXIT_FAILURE);
    }

    char str[200];
    while (fscanf(file,"%s",str) != EOF)
    {
        if (strcmp(str,"POINTS") == 0) break;
    }
    uint32_t num_points;
    fscanf(file,"%u %s",&num_points,str);
    for (uint32_t i = 0; i < num_points; i++)
    {
        struct network_point p;
        p.id = i;
        fscanf(file,"%lf %lf %lf",&p.x,&p.y,&p.z);

        points.push_back(p);
    }
    while (fscanf(file,"%s",str) != EOF)
    {
        if (strcmp(str,"LINES") == 0) break;
    }
    uint32_t num_lines, trash;
    fscanf(file,"%u %u",&num_lines,&trash);
    for (uint32_t i = 0; i < num_lines; i++)
    {
        struct network_line l;
        fscanf(file,"%u %u %u",&trash,&l.src,&l.dest);

        lines.push_back(l);
    }

    fclose(file);
}

void calculate_unitary_vector(const double x1, const double y1, const double z1, const double x2, const double y2, const double z2,\
                                double u[], double &norm)
{
    norm = calculate_euclidean_norm(x1,y1,z1,x2,y2,z2);
    u[0] = (x2-x1)/norm;
    u[1] = (y2-y1)/norm;
    u[2] = (z2-z1)/norm;   
}