#include "user_options.h"

struct user_options* new_user_options (int argc, char *argv[])
{
    struct user_options *result = (struct user_options*)malloc(sizeof(struct user_options));

    result->use_initial_network = false;

    read_config_file(result,argv[1]);

    return result;    
}

void free_user_options (struct user_options *the_options)
{
    if (the_options->walker_config)
        free_walker_config(the_options->walker_config);
    
    free(the_options);
}

void read_config_file (struct user_options *the_options, const char filename[])
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
    if(ini_parse(filename, parse_config_file, the_options) < 0) 
    {
        fprintf(stderr, "Error: Can't load the config file %s\n", filename);
        exit(EXIT_FAILURE);
    }

    // Set the function pointer to move and generate Walkers 
    set_walker_functions(the_options->walker_config);

    printf("%s\n",PRINT_LINE);

    fclose(file);
}

int parse_config_file(void *user, const char *section, const char *name, const char *value) 
{
    struct user_options *pconfig = (struct user_options*)user;

    if (SECTION_STARTS_WITH(MAIN_SECTION))
    {
        if (MATCH_NAME("number_of_iterations"))
        {
            pconfig->max_num_iter = (int)strtol(value, NULL, 10);
        }
        else if (MATCH_NAME("number_of_walker"))
        {
            pconfig->max_num_walkers = (int)strtol(value, NULL, 10);
        }
        else if (MATCH_NAME("seed"))
        {
            pconfig->seed = (uint32_t)strtol(value, NULL, 10);
        }
        else if (MATCH_NAME("root_pos_x"))
        {
            pconfig->root_pos[0] = (double)strtod(value, NULL);
        }
        else if (MATCH_NAME("root_pos_y"))
        {
            pconfig->root_pos[1] = (double)strtod(value, NULL);
        }
        else if (MATCH_NAME("root_pos_z"))
        {
            pconfig->root_pos[2] = (double)strtod(value, NULL);
        }
        else if (MATCH_NAME("initial_network_filename"))
        {
            pconfig->use_initial_network = true;
            pconfig->initial_network_filename = strdup(value);
        }
    }
    else if (SECTION_STARTS_WITH(WALKER_SECTION))
    {
        if (!pconfig->walker_config)
        {
            pconfig->walker_config = new_walker_config();
        }

        if (MATCH_NAME("library_name"))
        {
            pconfig->walker_config->library_name = strdup(value);
        }
        else if (MATCH_NAME("mesh_filename"))
        {
            pconfig->walker_config->mesh_filename = strdup(value);
        }
        else if (MATCH_NAME("map_filename"))
        {
            pconfig->walker_config->map_filename = strdup(value);
        }
        else
        {
            std::string key(name);
        
            pconfig->walker_config->params->insert(std::pair<std::string,double>(key,atof(value)));
        }
    }

    return 1;
}

void print_user_options (struct user_options *the_options)
{
    printf("%s\n",PRINT_LINE);
    printf("[user_options] number_of_iterations = %u\n",the_options->max_num_iter);
    printf("[user_options] number_of_walker = %u\n",the_options->max_num_walkers);
    printf("%s\n",PRINT_DOTS);
    if (the_options->walker_config)
        print_walker_config(the_options->walker_config);
    printf("%s\n",PRINT_LINE);
}