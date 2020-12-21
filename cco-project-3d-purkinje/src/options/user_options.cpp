#include "user_options.h"

User_Options::User_Options (const char filename[])
{
    this->use_cloud_points = false;
    this->use_local_optimization = false;
    this->use_obstacle = false;
    this->use_only_murray = false;
    this->use_pruning = false;
    this->use_pmj_location = false;
    this->use_initial_network = false;
    this->start_radius = -1;
    this->seed = 1;                             // Default value
    this->max_rand_offset = 1;                  // Default value
    this->gamma = 3.0;                          // Default value

    this->cost_function_config = NULL;
    this->local_opt_config = NULL;
    this->pmj_config = NULL;

    read_config_file(filename);
}

User_Options::~User_Options ()
{
    if (this->cost_function_config)
        delete this->cost_function_config;
    if (this->local_opt_config)
        delete this->local_opt_config;
}

void User_Options::read_config_file (const char filename[])
{
    printf("%s\n",PRINT_LINE);
    printf("[user_options] Reading configuration file:> \"%s\"\n",filename);

    // Open the config file for reading
    FILE *file = fopen(filename,"r");
    if (!file)
    {
        fprintf(stderr,"[-] Error reading configuration file '%s'\n",filename);
        exit(EXIT_FAILURE);
    }
    
    // Here we parse the config file
    if(ini_parse(filename, parse_config_file, this) < 0) 
    {
        fprintf(stderr, "Error: Can't load the config file %s\n", filename);
        exit(EXIT_FAILURE);
    }

    // Set the pruning function pointer
    //if (options->pruning_config)
    //    set_pruning_function(options->pruning_config);

    printf("%s\n",PRINT_LINE);

    fclose(file);

    // DEBUG
    //print();
}

int parse_config_file (void *user, const char *section, const char *name, const char *value)
{
    User_Options *pconfig = (User_Options*)user;

    if (SECTION_STARTS_WITH(MAIN_SECTION))
    {
        if (MATCH_NAME("N_term"))
        {
            pconfig->n_term = (int)strtol(value, NULL, 10);
        }
        else if (MATCH_NAME("root_x"))
        {
            pconfig->root_pos[0] = strtof(value, NULL);
        }
        else if (MATCH_NAME("root_y"))
        {
            pconfig->root_pos[1] = strtof(value, NULL);
        }
        else if (MATCH_NAME("root_z"))
        {
            pconfig->root_pos[2] = strtof(value, NULL);
        }
        else if (MATCH_NAME("seed"))
        {
            pconfig->seed = (uint32_t)strtol(value, NULL, 10);
        }
        else if (MATCH_NAME("max_rand_offset"))
        {
            pconfig->max_rand_offset = (uint32_t)strtol(value, NULL, 10);
        }
        else if (MATCH_NAME("gamma"))
        {
            pconfig->gamma = strtof(value, NULL);
        }
        else if (MATCH_NAME("lat_offset"))
        {
            pconfig->lat_offset = strtof(value, NULL);
        }
        else if (MATCH_NAME("use_only_murray"))
        {
            if (strcmp(value,"true") == 0 || strcmp(value,"yes") == 0)
                pconfig->use_only_murray = true;
            else if (strcmp(value,"false") == 0 || strcmp(value,"no") == 0)
                pconfig->use_only_murray = false;
            else
            {
                fprintf(stderr,"[user_options] Error reading configuration file! Invalid option in \"main\" section\n");
                exit(EXIT_FAILURE);
            }
        }
        else if (MATCH_NAME("start_radius"))
        {
            pconfig->start_radius = strtof(value, NULL);
        }
        else if (MATCH_NAME("use_initial_network"))
        {
            if (strcmp(value,"true") == 0 || strcmp(value,"yes") == 0)
                pconfig->use_initial_network = true;
            else if (strcmp(value,"false") == 0 || strcmp(value,"no") == 0)
                pconfig->use_initial_network = false;
            else
            {
                fprintf(stderr,"[user_options] Error reading configuration file! Invalid option in \"main\" section\n");
                exit(EXIT_FAILURE);
            }
        }
        else if (MATCH_NAME("initial_network_filename"))
        {
            pconfig->initial_network_filename = value;
        }
    }
    else if (SECTION_STARTS_WITH(SAVE_NETWORK_SECTION))
    {
        if (MATCH_NAME("output_dir"))
        {
            pconfig->output_dir = value;
        }
    }
    else if (SECTION_STARTS_WITH(CLOUD_SECTION))
    {
        if (MATCH_NAME("use_cloud_points"))
        {
            if (strcmp(value,"true") == 0 || strcmp(value,"yes") == 0)
                pconfig->use_cloud_points = true;
            else if (strcmp(value,"false") == 0 || strcmp(value,"no") == 0)
                pconfig->use_cloud_points = false;
            else
            {
                fprintf(stderr,"[user_options] Error reading configuration file! Invalid option in \"cloud_points\" section\n");
                exit(EXIT_FAILURE);
            }
        }
        else if (MATCH_NAME("cloud_points_filename"))
        {
            pconfig->cloud_points_filename = value;
        }
        else if (MATCH_NAME("use_obstacle"))
        {
            if (strcmp(value,"true") == 0 || strcmp(value,"yes") == 0)
                pconfig->use_obstacle = true;
            else if (strcmp(value,"false") == 0 || strcmp(value,"no") == 0)
                pconfig->use_obstacle = false;
            else
            {
                fprintf(stderr,"[user_options] Error reading configuration file! Invalid option in \"cloud_points\" section\n");
                exit(EXIT_FAILURE);
            }
        }
        else if (MATCH_NAME("obstacle_filename"))
        {
            pconfig->obstacle_filename = value;
        }
    }
    else if (SECTION_STARTS_WITH(PMJ_SECTION))
    {
        if (!pconfig->pmj_config)
        {
            pconfig->pmj_config = new PMJConfig();
        }

        if (MATCH_NAME("use_pmj_location"))
        {
            if (strcmp(value,"true") == 0 || strcmp(value,"yes") == 0)
                pconfig->use_pmj_location = true;
            else if (strcmp(value,"false") == 0 || strcmp(value,"no") == 0)
                pconfig->use_pmj_location = false;
            else
            {
                fprintf(stderr,"[user_options] Error reading configuration file! Invalid option in \"cloud_points\" section\n");
                exit(EXIT_FAILURE);
            }
        }
        else if (MATCH_NAME("pmj_location_filename"))
        {
            pconfig->pmj_config->location_filename = value;
        }
        else if (MATCH_NAME("max_pmj_connection_tries"))
        {
            pconfig->pmj_config->max_connection_tries = (uint32_t)strtol(value, NULL, 10);
        }
        else if (MATCH_NAME("pmj_connection_rate"))
        {
            pconfig->pmj_config->connection_rate = (uint32_t)strtol(value, NULL, 10);
        }
        else if (MATCH_NAME("pmj_region_radius"))
        {
            pconfig->pmj_config->region_radius = strtof(value, NULL);
        }
        else if (MATCH_NAME("lat_error_tolerance"))
        {
            pconfig->pmj_config->lat_error_tolerance = strtof(value, NULL);
        }
    }
    else if (SECTION_STARTS_WITH(LOCAL_OPT_SECTION))
    {
        if (!pconfig->local_opt_config)
        {
            pconfig->local_opt_config = new LocalOptimizationConfig();
        }

        if (MATCH_NAME("use_local_optimization"))
        {
            if (strcmp(value,"true") == 0 || strcmp(value,"yes") == 0)
                pconfig->use_local_optimization = true;
            else
                pconfig->use_local_optimization = false;
        }
        else if (MATCH_NAME("local_optimization_function"))
        {
            pconfig->local_opt_config->function_name = value;
        }
        else if (MATCH_NAME("library_name"))
        {
            pconfig->local_opt_config->library_name = value;
        }
    }
    else if (SECTION_STARTS_WITH(COST_FUNCTION_SECTION))
    {
        if (!pconfig->cost_function_config)
        {
            pconfig->cost_function_config = new CostFunctionConfig();
        }

        if (MATCH_NAME("library_name"))
        {
            pconfig->cost_function_config->library_name = value;
        }   
        else if (MATCH_NAME("function_name"))
        {
            pconfig->cost_function_config->function_name = value;
        }
        else
        {
            std::string key(name);
        
            pconfig->cost_function_config->params->insert(std::pair<std::string,double>(key,atof(value)));
        }
    }
    /*
    else if (SECTION_STARTS_WITH(PRUNING_SECTION))
    {
        if (!pconfig->pruning_config)
        {
            pconfig->pruning_config = new_pruning_config();
        }

        if (MATCH_NAME("use_pruning"))
        {
            if (strcmp(value,"true") == 0 || strcmp(value,"yes") == 0)
                pconfig->use_pruning = true;
            else
                pconfig->use_pruning = false;
        }
        else if (MATCH_NAME("pruning_function"))
        {
            pconfig->pruning_config->function_name = strdup(value);
        }
        else if (MATCH_NAME("library_name"))
        {
            pconfig->pruning_config->library_name = strdup(value);
        }
        else
        {
            std::string key(name);
        
            pconfig->pruning_config->params->insert(std::pair<std::string,double>(key,atof(value)));
        }
    }
    */

    return 1;

}

void User_Options::print ()
{
    printf("********************* user_options *********************\n");
    printf("N_term = %u\n",this->n_term);
    printf("seed = %u\n",this->seed);
    printf("max_rand_offset = %u\n",this->max_rand_offset);
    printf("root_pos[3] = { %g , %g , %g }\n",this->root_pos[0],this->root_pos[1],this->root_pos[2]);
    printf("gamma = %g\n",this->gamma);
    printf("start_radius = %g mm\n",this->start_radius);
    printf("output_dir = %s\n",this->output_dir.c_str());

    printf("%s\n",PRINT_DOTS);
    if (this->use_cloud_points)
        printf("cloud_points_filename = %s\n",this->cloud_points_filename.c_str());
    else
        printf("cloud_points_filename = NULL\n");
    if (this->use_obstacle)
        printf("obstacle_filename = %s\n",this->obstacle_filename.c_str());
    else
        printf("obstacle_filename = NULL\n");
    if (this->use_initial_network)
        printf("Initial network filename = %s\n",this->initial_network_filename.c_str());
    else
        printf("Initial network filename = NULL\n");
    printf("%s\n",PRINT_DOTS);
    this->cost_function_config->print();
    printf("%s\n",PRINT_DOTS);
    if (this->use_local_optimization)
        this->local_opt_config->print();
    else
        printf("Local optimization = FALSE\n");
    if (this->use_pmj_location)
        this->pmj_config->print();
    else
        printf("PMJ location = FALSE\n");
    printf("********************************************************\n");
}