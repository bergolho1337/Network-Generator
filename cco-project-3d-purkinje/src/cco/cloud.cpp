#include "cloud.h"

Cloud::Cloud ()
{
    this->filename = "";
    this->cur_index = 0;
}

Cloud::Cloud (CloudConfig *config)
{
    this->cur_index = 0;
    this->filename = config->filename;
    bool sucess = read_points_from_vtk(config->filename.c_str(),this->points);
    if (sucess) 
    {
        printf("[cloud] Cloud of points was sucessfully loaded from: '%s'!\n",config->filename.c_str());
    }
    else
    {
        fprintf(stderr,"[-] ERROR! Cannot read cloud of points file: '%s'!\n",config->filename.c_str());
        exit(EXIT_FAILURE);
    }
    this->connected.assign(this->points.size(),false);
    this->aprox.assign(this->points.size(),__DBL_MAX__);
    this->error.assign(this->points.size(),__DBL_MIN__);
    this->ref.assign(this->points.size(),0);
    for (uint32_t i = 0; i < this->points.size(); i++) this->ref[i] = this->points[i]->lat;
    
    //print();
}

uint32_t Cloud::sort_point (Point *p, const uint32_t rand_offset)
{
    uint32_t offset = rand() % rand_offset + 1;

    // Reset the counter
    if (this->cur_index > this->points.size()-1) this->cur_index = this->cur_index % (this->points.size()-1);

    uint32_t selected_index = this->cur_index;
    
    // Get the current cloud point data
    uint32_t id;
    bool is_active;
    double pos[3], ref_lat;
    id = selected_index;
    pos[0] = this->points[selected_index]->x;
    pos[1] = this->points[selected_index]->y;
    pos[2] = this->points[selected_index]->z;
    ref_lat = this->points[selected_index]->lat;
    is_active = this->points[selected_index]->is_active;

    // Set the point data
    p->setId(id);
    p->setCoordinate(pos);
    p->setLAT(ref_lat);
    p->setActive(is_active);
    
    // Increase the counter for the next cloud point
    this->cur_index += offset;

    return selected_index;
}

Cloud::~Cloud ()
{
    for (uint32_t i = 0; i < this->points.size(); i++)
        delete this->points[i];
}

Cloud* Cloud::copy ()
{
    Cloud *result = new Cloud();
    result->filename = this->filename;
    for (uint32_t i = 0; i < this->points.size(); i++)
    {
        result->points.push_back(this->points[i]);
        result->connected.push_back(this->connected[i]);
        result->error.push_back(this->error[i]);
        result->aprox.push_back(this->aprox[i]);
        result->ref.push_back(this->ref[i]);
    }
    return result;
}

Cloud* Cloud::get_points_around_region (Point *pmj_point, const uint32_t num_points, const double region_radius)
{
    Cloud *result = new Cloud();
    double r = region_radius;
    while (result->points.size() < num_points)
    {
        for (uint32_t i = 0; i < this->points.size(); i++)
        {
            double dist = euclidean_norm(this->points[i]->x,this->points[i]->y,this->points[i]->z,\
                                        pmj_point->x,pmj_point->y,pmj_point->z);
            if (dist < r)
            {
                Point *p = this->points[i]->copy();
                result->points.push_back(p);
            }
        }
        r *= 1.2;
    }
    return result;
}

void Cloud::concatenate (Cloud *input)
{
    uint32_t offset_points = this->points.size();
    for (uint32_t i = 0; i < input->points.size(); i++)
    {
        Point *p = new Point(input->points[i]);
        bool connected = input->connected[i];
        double error = input->error[i];
        double aprox = input->aprox[i];
        double ref = input->ref[i];

        p->id += offset_points;
        
        this->points.push_back(p);
        this->connected.push_back(connected);
        this->error.push_back(error);
        this->aprox.push_back(aprox);
        this->ref.push_back(ref);
    }
}

void Cloud::print ()
{
    printf("[cloud] Cloud points filename: '%s'\n",this->filename.c_str());
    printf("[cloud] Number of points: %u\n",this->points.size());
    for (uint32_t i = 0; i < this->points.size(); i++)
        this->points[i]->print();
    for (uint32_t i = 0; i < this->points.size(); i++)
        printf("Point %u -- [Ref: %g] [Aprox: %g] [Error:%g]\n",i,ref[i],aprox[i],error[i]);
}