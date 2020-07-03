#include "cloud.h"

Cloud_Point::Cloud_Point (const char filename[])
{
    FILE *file = fopen(filename,"r");

    uint32_t num_points;
    fscanf(file,"%u",&num_points);

    for (uint32_t i = 0; i < num_points; i++)
    {
        double pos[3];
        fscanf(file,"%lf %lf %lf",&pos[0],&pos[1],&pos[2]);

        Point_3d p(i,pos[0],pos[1],pos[2]);
        this->the_points.push_back(p);
        this->taken.push_back(false);
    }

    fclose(file);
}

void Cloud_Point::print ()
{
    for (uint32_t i = 0; i < this->the_points.size(); i++)
    {
        printf("Point %u -- (%g,%g,%g)\n",this->the_points[i].id,this->the_points[i].x,this->the_points[i].y,this->the_points[i].z);
    }
    printf("Number of points = %u\n",this->the_points.size());
}