#include "pmj.h"

PMJ::PMJ ()
{
    this->total_num_connected = 0;
    this->max_connection_tries = 100;
    this->connection_rate = __UINT32_MAX__;
    this->region_radius = 0.002;
    this->lat_error_tolerance = 2.0;
    this->location_filename = "";
}

PMJ::PMJ (PMJConfig *config)
{
    this->total_num_connected = 0;
    this->max_connection_tries = config->max_connection_tries;
    this->connection_rate = config->connection_rate;
    this->region_radius = config->region_radius;
    this->lat_error_tolerance = config->lat_error_tolerance;
    this->location_filename = config->location_filename;

    bool sucess = read_points_from_vtk(this->location_filename.c_str(),this->points);
    if (sucess) 
    {
        printf("[pmj] PMJ location file was sucessfully loaded from: '%s'!\n",this->location_filename.c_str());
    }
    else
    {
        fprintf(stderr,"[-] ERROR! Cannot read PMJ location file: '%s'!\n",config->location_filename.c_str());
        exit(EXIT_FAILURE);
    }

    uint32_t n = this->points.size();
    this->connected.assign(n,false);
    this->error.assign(n,__DBL_MAX__);
    this->aprox.assign(n,0);

    // DEBUG
    //print();
}

PMJ::~PMJ ()
{
    for (uint32_t i = 0; i < this->points.size(); i++)
        delete this->points[i];
}

PMJ* PMJ::copy ()
{
    PMJ *result = new PMJ();
    result->total_num_connected = this->total_num_connected;
    result->max_connection_tries = this->max_connection_tries;
    result->connection_rate = this->connection_rate;
    result->region_radius = this->region_radius;
    result->lat_error_tolerance = this->lat_error_tolerance;
    result->location_filename = this->location_filename;
    for (uint32_t i = 0; i < this->points.size(); i++)
    {
        result->points.push_back(this->points[i]);
        result->connected.push_back(this->connected[i]);
        result->aprox.push_back(this->aprox[i]);
        result->error.push_back(this->error[i]);
    }   
    return result;
}

void PMJ::print ()
{
    printf("[pmj] Total number of PMJ's to connect = %u\n",this->points.size());
    printf("[pmj] Maximum number of connection tries = %u\n",this->max_connection_tries);
    printf("[pmj] Connection rate = %u\n",this->connection_rate);
    printf("[pmj] Region radius = %g m\n",this->region_radius);
    printf("[pmj] LAT error tolerance = %g ms\n",this->lat_error_tolerance);
    printf("[pmj] PMJ location filename = %s\n",this->location_filename.c_str());
    for (uint32_t i = 0; i < this->points.size(); i++)
        this->points[i]->print();
}
