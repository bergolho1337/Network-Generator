#include "pmj.h"

PMJ::PMJ ()
{
    initialize_parameters();
}

PMJ::PMJ (PMJConfig *config)
{
    if (config->using_pmj)
    {
        initialize_parameters(config);
        initialize_points(config);
        initialize_arrays(config);
        //initialize_packages();
        initialize_packages_2();
    }
    else
    {
        printf("[pmj] No PMJ configuration provided\n");
        initialize_parameters();
    }
    
    // DEBUG
    //print();
}

PMJ::~PMJ ()
{
    if (this->cost_fn)
        delete this->cost_fn;
    for (uint32_t i = 0; i < this->points.size(); i++)
        delete this->points[i];
}

void PMJ::initialize_parameters ()
{
    this->total_num_connected = 0;
    this->package_size = 10;
    this->package_head = 0;
    this->cur_package = 0;
    this->max_connection_tries = 100;
    this->connection_rate = __UINT32_MAX__;
    this->region_radius = 0.002;
    this->lat_error_tolerance = 2.0;
    this->location_filename = "";
    this->cost_fn = new MaximizeCustomFunction;
}

void PMJ::initialize_parameters (PMJConfig *config)
{
    this->total_num_connected = 0;
    this->package_size = 10;
    this->package_head = 0;
    this->cur_package = 0;
    this->max_connection_tries = config->max_connection_tries;
    this->connection_rate = config->connection_rate;
    this->region_radius = config->region_radius;
    this->lat_error_tolerance = config->lat_error_tolerance;
    this->location_filename = config->location_filename;
    this->cost_fn = new MaximizeCustomFunction;
}

void PMJ::initialize_arrays (PMJConfig *config)
{
    uint32_t n = this->points.size();
    for (uint32_t i = 0; i < n; i++) this->reference.push_back(this->points[i]->lat);
    this->connected.assign(n,false);
    this->error.assign(n,__DBL_MAX__);
    this->aprox.assign(n,0);
    this->penalty.assign(n,0);
}

void PMJ::initialize_points (PMJConfig *config)
{
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
}

void PMJ::initialize_packages ()
{
    this->packages.assign(this->package_size,std::vector<Point*>());
    std::sort(this->points.begin(),this->points.end(),comparePoint);

    double offset = (this->points.back()->lat-this->points.front()->lat) / (double)this->package_size;
    double base = this->points.front()->lat;
    for (uint32_t i = 0; i < this->package_size; i++)
    {
        double min_value = base + i*offset;
        double max_value = min_value + offset + 1.0e-8;

        uint32_t counter = 0;
        for (uint32_t j = 0; j < this->points.size(); j++)
        {
            double value = this->points[j]->lat;
            if (value >= min_value && value <= max_value)
                this->packages[i].push_back(this->points[j]);
        }
    }
}

void PMJ::initialize_packages_2 ()
{
    std::sort(this->points.begin(),this->points.end(),comparePoint);

    if (this->points.size() < this->package_size) this->package_size = this->points.size();
    for (uint32_t i = 0; i < this->package_size; i++)
    {
        this->package.push_back(i);
        this->package_head++;
    }
}

void PMJ::concatenate (PMJ *input)
{
    this->total_num_connected += input->total_num_connected;
        
    uint32_t offset_points = this->points.size();
    for (uint32_t i = 0; i < input->points.size(); i++)
    {
        Point *p = new Point(input->points[i]);
        bool connected = input->connected[i];
        double aprox = input->aprox[i];
        double reference = input->reference[i];
        double error = input->error[i];

        p->id += offset_points;

        this->points.push_back(p);
        this->connected.push_back(connected);
        this->aprox.push_back(aprox);
        this->reference.push_back(reference);
        this->error.push_back(error);
    }
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
    result->package_size = this->package_size;
    result->package_head = 0;
    result->cost_fn = new MaximizeCustomFunction;
    
    for (uint32_t i = 0; i < this->points.size(); i++)
    {
        result->points.push_back(this->points[i]);
        result->connected.push_back(this->connected[i]);
        result->reference.push_back(this->reference[i]);
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

bool comparePoint (Point *a, Point *b)
{
    return a->lat < b->lat;
}