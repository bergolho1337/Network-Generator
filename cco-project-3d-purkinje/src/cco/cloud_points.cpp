#include "cloud_points.h"

Cloud_Point::Cloud_Point (std::string filename)
{
    bool sucess = read_points_from_vtk(this->filename.c_str(),this->points);
    if (sucess) 
    {
        printf("[cco] Cloud of points was sucessfully loaded from: '%s'!\n",this->filename.c_str());
    }
    else
    {
        fprintf(stderr,"[-] ERROR! Cannot read cloud of points file: '%s'!\n",filename.c_str());
        exit(EXIT_FAILURE);
    }
}

Cloud_Point::~Cloud_Point ()
{
    for (uint32_t i = 0; i < this->points.size(); i++)
        delete this->points[i];
}