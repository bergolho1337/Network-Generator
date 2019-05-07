#include "user_options.h"

struct user_options* new_user_options (int argc, char *argv[])
{
    struct user_options *result = (struct user_options*)malloc(sizeof(struct user_options));

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

    FILE *file = fopen(filename,"r");
    
    // Reading [main] section
    read_main_section(options,file);

    // Reading [cloud_points] section
    read_cloud_points_section(options,file);

    // Reading [local_optimization] section
    read_local_optimization_section(options,file);

    // Reading [cost_function] section
    read_cost_function_section(options,file);

    printf("%s\n",PRINT_LINE);

    fclose(file);
}

void read_main_section (struct user_options *options, FILE *file)
{
    char str[MAX_FILENAME_SIZE] = "", trash[MAX_FILENAME_SIZE] = "", value[MAX_FILENAME_SIZE] = "";

    while (strcmp(str,"[main]") != 0)
        fscanf(file,"%s",str);

    fscanf(file,"%s %s %lf",str,trash,&options->Q_perf);
    fscanf(file,"%s %s %lf",str,trash,&options->p_perf);
    fscanf(file,"%s %s %lf",str,trash,&options->p_term);
    fscanf(file,"%s %s %lf",str,trash,&options->r_perf);
    fscanf(file,"%s %s %d",str,trash,&options->N_term);
}

void read_cloud_points_section (struct user_options *options, FILE *file)
{
    char str[MAX_FILENAME_SIZE] = "", trash[MAX_FILENAME_SIZE] = "", value[MAX_FILENAME_SIZE] = "";

    while (strcmp(str,"[cloud_points]") != 0)
        fscanf(file,"%s",str);

    fscanf(file,"%s %s %s",str,trash,value);
    if (strcmp(value,"true") == 0 || strcmp(value,"yes") == 0)
    {
        options->use_cloud_points = true;
        fscanf(file,"%s %s %s",str,trash,options->cloud_points_filename);
    }
    else if (strcmp(value,"false") == 0 || strcmp(value,"no") == 0)
    {
        options->use_cloud_points = false;
    }
    else
    {
        fprintf(stderr,"[user_options] Error reading configuration file! Invalid option in \"cloud_points\" section\n");
        exit(EXIT_FAILURE);
    }
}

void read_cost_function_section (struct user_options *options, FILE *file)
{
    char str[MAX_FILENAME_SIZE] = "", trash[MAX_FILENAME_SIZE] = "", value[MAX_FILENAME_SIZE] = "";

    while (strcmp(str,"[cost_function]") != 0)
        fscanf(file,"%s",str);
    
    options->config = new_cost_function_config();

    // Read cost function library name
    fscanf(file,"%s %s %s",str,trash,value);
    options->config->library_name = (char*)malloc(sizeof(char)*(strlen(value)+1));
    strcpy(options->config->library_name,value);

    // Read cost function name
    fscanf(file,"%s %s %s",str,trash,value);
    options->config->function_name = (char*)malloc(sizeof(char)*(strlen(value)+1));
    strcpy(options->config->function_name,value);

    // Read cost function parameters
    while (fscanf(file,"%s %s %s",str,trash,value) != EOF)
    {   

        std::string key(str);
        
        options->config->params->insert(std::pair<std::string,double>(key,atof(value)));
    }

    // Set the cost function pointer
    set_cost_function(options->config);

    //options->config->function(NULL,options->config);

}

void read_local_optimization_section (struct user_options *options, FILE *file)
{
    char str[MAX_FILENAME_SIZE] = "", trash[MAX_FILENAME_SIZE] = "", value[MAX_FILENAME_SIZE] = "";

    while (strcmp(str,"[local_optimization]") != 0)
        fscanf(file,"%s",str);
    
    options->local_opt_config = new_local_optimization_config();

    // Read if we are going to use or not local optimization
    fscanf(file,"%s %s %s",str,trash,value);
    if (strcmp(value,"true") == 0 || strcmp(value,"yes") == 0)
    {
        options->use_local_optimization = true;

        // The next input should be the name of the local optimization function
        fscanf(file,"%s %s %s",str,trash,value);
        options->local_opt_config->name = (char*)malloc(sizeof(char)*(strlen(value)+1));
        strcpy(options->local_opt_config->name,value);

        // Set the local optimization function pointer
        set_local_optimization_function(options->local_opt_config);

        // Call the function for testing
        //std::vector<struct point*> test_pos;
        //options->local_opt_config->function(NULL,NULL,NULL,test_pos);
    }
    else if (strcmp(value,"false") == 0 || strcmp(value,"no") == 0)
    {
        options->use_local_optimization = false;
    }
    else
    {
        fprintf(stderr,"[local_optimization] Error reading configuration file! Invalid option in \"local_optimization\" section\n");
        exit(EXIT_FAILURE);
    }

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