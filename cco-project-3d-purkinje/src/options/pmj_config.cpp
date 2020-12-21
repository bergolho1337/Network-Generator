#include "pmj_config.h"

PMJConfig::PMJConfig ()
{
    this->lat_error_tolerance = 2.0;
    this->max_connection_tries = 100;
    this->connection_rate = __UINT32_MAX__;
    this->region_radius = 0.002;
    this->location_filename = "";
}

PMJConfig::~PMJConfig ()
{

}

void PMJConfig::print ()
{
    printf("LAT tolerance error = %g ms\n",this->lat_error_tolerance);
    printf("Max. PMJ connection tries = %u\n",this->max_connection_tries);
    printf("PMJ connection rate = %u\n",this->connection_rate);
    printf("PMJ region radius = %g m\n",this->region_radius);
    printf("PMJ location filename = \"%s\"\n",this->location_filename.c_str());
}