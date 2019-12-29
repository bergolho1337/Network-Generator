#include "user_options.h"

struct user_options* new_user_options (int argc, char *argv[])
{
    struct user_options *result = (struct user_options*)malloc(sizeof(struct user_options));

    result->use_cloud_points = false;
    result->use_local_optimization = false;

    read_config_file(result,argv[1]);
    
    return result;
}

void free_user_options (struct user_options *options)
{
    free_cost_function_config(options->config);

    free_local_optimization_config(options->local_opt_config);
    
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
            pconfig->N_term = (int)strtol(value, NULL, 10);
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
            pconfig->Q_perf = strtof(value, NULL);
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
            pconfig->V_perf = strtof(value, NULL);
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

    return 1;
}

void print_user_options (struct user_options *options)
{
    printf("********************* user_options *********************\n");
    printf("Q_perf = %g\n",options->Q_perf);
    printf("p_perf = %g\n",options->p_perf);
    printf("p_term = %g\n",options->p_term);
    printf("r_perf = %g\n",options->r_perf);
    printf("N_term = %d\n",options->N_term);
    printf("%s\n",PRINT_DOTS);
    if (options->use_cloud_points)
        printf("cloud_points_filename = %s\n",options->cloud_points_filename);
    else
        printf("cloud_points_filename = NULL\n");
    printf("%s\n",PRINT_DOTS);
    print_cost_function_config(options->config);
    printf("%s\n",PRINT_DOTS);
    if (options->use_local_optimization)
        print_local_optimization_function_config(options->local_opt_config);
    else
        printf("Local optimization = FALSE\n");
    printf("********************************************************\n");
}