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
    this->print();
}

int parse_config_file(void *user, const char *section, const char *name, const char *value) 
{
    User_Options *pconfig = (User_Options*)user;

    if (SECTION_STARTS_WITH(MAIN_SECTION))
    {
        if (MATCH_NAME("max_iterations"))
        {
            pconfig->max_iterations = (uint32_t)strtol(value, NULL, 10);
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
        else if (MATCH_NAME("initial_length"))
        {
            pconfig->initial_length = strtof(value, NULL);
        }
        else if (MATCH_NAME("initial_angle"))
        {
            pconfig->initial_angle = strtof(value, NULL);
        }
        else if (MATCH_NAME("initial_diameter"))
        {
            pconfig->initial_diameter = strtof(value, NULL);
        }
        else if (MATCH_NAME("length_decrease_ratio"))
        {
            pconfig->length_decrease_ratio = strtof(value, NULL);
        }
        else if (MATCH_NAME("angle_decrease_ratio"))
        {
            pconfig->angle_decrease_ratio = strtof(value, NULL);
        }
        else if (MATCH_NAME("diameter_decrease_ratio"))
        {
            pconfig->diameter_decrease_ratio = strtof(value, NULL);
        }
    }

    return 1;
}

void User_Options::print ()
{
    printf("********************* user_options *********************\n");
    printf("max_iterations = %u\n",this->max_iterations);
    printf("initial_length = %g\n",this->initial_length);
    printf("initial_angle = %g\n",this->initial_angle);
    printf("length_decrease_ratio = %g\n",this->length_decrease_ratio);
    printf("angle_decrease_ratio = %g\n",this->angle_decrease_ratio);
    printf("root_pos = (%g,%g,%g)\n",this->root_pos[0],this->root_pos[1],this->root_pos[2]);
    printf("********************************************************\n");
}