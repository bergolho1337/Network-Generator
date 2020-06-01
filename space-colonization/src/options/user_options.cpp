#include "user_options.h"

User_Options::User_Options (int argc, char *argv[])
{
    read_config_file(argv[1]);
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

    printf("%s\n",PRINT_LINE);

    fclose(file);

    // DEBUG
    //print();
}

int parse_config_file(void *user, const char *section, const char *name, const char *value) 
{
    User_Options *pconfig = (User_Options*)user;

    if (SECTION_STARTS_WITH(MAIN_SECTION))
    {
        if (MATCH_NAME("number_of_iterations"))
        {
            pconfig->number_of_iterations = (uint32_t)strtol(value, NULL, 10);
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
        else if (MATCH_NAME("root_dir_x"))
        {
            pconfig->root_dir_pos[0] = strtof(value, NULL);
        }
        else if (MATCH_NAME("root_dir_y"))
        {
            pconfig->root_dir_pos[1] = strtof(value, NULL);
        }
        else if (MATCH_NAME("root_dir_z"))
        {
            pconfig->root_dir_pos[2] = strtof(value, NULL);
        }
        else if (MATCH_NAME("segment_length"))
        {
            pconfig->segment_length = strtof(value, NULL);
        }
        else if (MATCH_NAME("leaves_offset"))
        {
            pconfig->leaves_offset = (uint32_t)strtol(value, NULL, 10);
        }
    }
    else if (SECTION_STARTS_WITH(SPACE_COLONIZATION_SECTION))
    {
        if (MATCH_NAME("min_distance"))
        {
            pconfig->min_distance = strtof(value, NULL);
        }
        else if (MATCH_NAME("max_distance"))
        {
            pconfig->max_distance = strtof(value, NULL);
        }
    }
    else if (SECTION_STARTS_WITH(CLOUD_SECTION))
    {
        if (MATCH_NAME("cloud_points_filename"))
        {
            pconfig->cloud_points_filename = strdup(value);
        }
    }

    return 1;
}

void User_Options::print ()
{
    printf("********************* user_options *********************\n");
    printf("number_of_iterations = %u\n",this->number_of_iterations);
    printf("root_pos = (%g,%g,%g)\n",this->root_pos[0],this->root_pos[1],this->root_pos[2]);
    printf("root_dir = (%g,%g,%g)\n",this->root_dir_pos[0],this->root_dir_pos[1],this->root_dir_pos[2]);
    printf("segment_length = %g\n",this->segment_length);
    printf("leaves_offset = %u\n",this->leaves_offset);
    printf("--------------------------------------------------------\n");
    printf("cloud_points_filename = %s\n",this->cloud_points_filename);
    printf("********************************************************\n");
}