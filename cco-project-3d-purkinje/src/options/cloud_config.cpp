#include "cloud_config.h"

CloudConfig::CloudConfig ()
{
    this->filename = "";
    this->rand_offset = 2;
}

CloudConfig::~CloudConfig ()
{

}

void CloudConfig::print ()
{
    printf("cloud_points_filename = %s\n",this->filename.c_str());
}