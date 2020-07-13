#include "user_options.h"

struct user_options* new_user_options (int argc, char *argv[])
{
    struct user_options *result = (struct user_options*)malloc(sizeof(struct user_options));

    result->use_cloud_points = false;
    result->use_local_optimization = false;
    result->use_obstacle = false;
    result->use_only_murray = false;
    result->use_pruning = false;
    result->use_pmj_location = false;
    result->use_initial_network = false;
    result->start_radius = -1;
    result->seed = 1;                           // Default value
    result->max_rand_offset = 1;                // Default value
    result->gamma = 3.0;                        // Default value

    read_config_file(result,argv[1]);
    
    return result;
}

void free_user_options (struct user_options *options)
{
    free_cost_function_config(options->config);

    free_local_optimization_config(options->local_opt_config);
    
    free(options->output_dir);
    free(options);
}

void read_config_file (struct user_options *options, const char filename[])
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
    if(ini_parse(filename, parse_config_file, options) < 0) 
    {
        fprintf(stderr, "Error: Can't load the config file %s\n", filename);
        exit(EXIT_FAILURE);
    }

    // Set the cost function pointer
    set_cost_function(options->config);

    // Set the local optimization function pointer
    if (options->use_local_optimization)
        set_local_optimization_function(options->local_opt_config);

    printf("%s\n",PRINT_LINE);

    fclose(file);

    //print_user_options(options);
}

int parse_config_file(void *user, const char *section, const char *name, const char *value) 
{
    struct user_options *pconfig = (struct user_options *)user;

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
        else if (MATCH_NAME("Q_perf"))
        {
            pconfig->q_perf = strtof(value, NULL);
        }
        else if (MATCH_NAME("p_perf"))
        {
            pconfig->p_perf = strtof(value, NULL);
        }
        else if (MATCH_NAME("p_term"))
        {
            pconfig->p_term = strtof(value, NULL);
        }
        else if (MATCH_NAME("V_perf"))
        {
            pconfig->v_perf = strtof(value, NULL);
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
            pconfig->initial_network_filename = strdup(value);
        }
    }
    else if (SECTION_STARTS_WITH(SAVE_NETWORK_SECTION))
    {
        if (MATCH_NAME("output_dir"))
        {
            pconfig->output_dir = strdup(value);
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
            pconfig->cloud_points_filename = strdup(value);
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
            pconfig->obstacle_filename = strdup(value);
        }
        else if (MATCH_NAME("use_pmj_location"))
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
            pconfig->pmj_location_filename = strdup(value);
        }
    }
    else if (SECTION_STARTS_WITH(LOCAL_OPT_SECTION))
    {
        if (!pconfig->local_opt_config)
        {
            pconfig->local_opt_config = new_local_optimization_config();
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
            pconfig->local_opt_config->name = strdup(value);
        }
    }
    else if (SECTION_STARTS_WITH(COST_FUNCTION_SECTION))
    {
        if (!pconfig->config)
        {
            pconfig->config = new_cost_function_config();
        }

        if (MATCH_NAME("library_name"))
        {
            pconfig->config->library_name = strdup(value);
        }   
        else if (MATCH_NAME("function_name"))
        {
            pconfig->config->function_name = strdup(value);
        }
        else
        {
            std::string key(name);
        
            pconfig->config->params->insert(std::pair<std::string,double>(key,atof(value)));
        }
    }
    else if (SECTION_STARTS_WITH(PRUNING_SECTION))
    {
        if (MATCH_NAME("use_pruning"))
        {
            if (strcmp(value,"true") == 0 || strcmp(value,"yes") == 0)
                pconfig->use_pruning = true;
            else
                pconfig->use_pruning = false;
        }
        else if (MATCH_NAME("A"))
        {
            pconfig->a = strtof(value, NULL);
        }
        else if (MATCH_NAME("B"))
        {
            pconfig->b = strtof(value, NULL);
        }
        else if (MATCH_NAME("C"))
        {
            pconfig->c = strtof(value, NULL);
        }
    }

    return 1;
}

void print_user_options (struct user_options *options)
{
    printf("********************* user_options *********************\n");
    printf("Q_perf = %g\n",options->q_perf);
    printf("p_perf = %g\n",options->p_perf);
    printf("p_term = %g\n",options->p_term);
    printf("V_perf = %g\n",options->v_perf);
    printf("N_term = %d\n",options->n_term);
    printf("%s\n",PRINT_DOTS);
    if (options->use_cloud_points)
        printf("cloud_points_filename = %s\n",options->cloud_points_filename);
    else
        printf("cloud_points_filename = NULL\n");
    if (options->use_obstacle)
        printf("obstacle_filename = %s\n",options->cloud_points_filename);
    else
        printf("obstacle_filename = NULL\n");
    printf("%s\n",PRINT_DOTS);
    print_cost_function_config(options->config);
    printf("%s\n",PRINT_DOTS);
    if (options->use_local_optimization)
        print_local_optimization_function_config(options->local_opt_config);
    else
        printf("Local optimization = FALSE\n");
    printf("********************************************************\n");
}