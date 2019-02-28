#include "user_options.h"

struct user_options* new_user_options (int argc, char *argv[])
{
    struct user_options *result = (struct user_options*)malloc(sizeof(struct user_options));

    read_config_file(result,argv[1]);
    
    return result;
}

void free_user_options (struct user_options *options)
{
    free(options);
}

void read_config_file (struct user_options *options, const char filename[])
{
    printf("%s\n",PRINT_LINE);
    printf("[user_options] Reading configuration file:> \"%s\"\n",filename);

    FILE *file = fopen(filename,"r");
    char str[MAX_FILENAME_SIZE] = "", trash[MAX_FILENAME_SIZE] = "", value[MAX_FILENAME_SIZE] = "";
    
    while (strcmp(str,"[main]") != 0)
        fscanf(file,"%s",str);

    // Reading [main] section
    fscanf(file,"%s %s %lf",str,trash,&options->Q_perf);
    fscanf(file,"%s %s %lf",str,trash,&options->p_perf);
    fscanf(file,"%s %s %lf",str,trash,&options->p_term);
    fscanf(file,"%s %s %lf",str,trash,&options->r_perf);
    fscanf(file,"%s %s %d",str,trash,&options->N_term);

    while (strcmp(str,"[cost_function]") != 0)
        fscanf(file,"%s",str);
    
    // Reading [cost_function] section
    fscanf(file,"%s %s %s",str,trash,options->cost_function_name);

    while (strcmp(str,"[cloud_points]") != 0)
        fscanf(file,"%s",str);

    // Reading [cloud_points] section
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
    printf("%s\n",PRINT_LINE);

    fclose(file);
}

void print_user_options (struct user_options *options)
{
    printf("********************* user_options *********************\n");
    printf("Q_perf = %g\n",options->Q_perf);
    printf("p_perf = %g\n",options->p_perf);
    printf("p_term = %g\n",options->p_term);
    printf("r_perf = %g\n",options->r_perf);
    printf("N_term = %d\n",options->N_term);
    printf("--------------------------------------------------------\n");
    printf("cost_function_name = %s\n",options->cost_function_name);
    printf("--------------------------------------------------------\n");
    if (options->use_cloud_points)
        printf("cloud_points_filename = %s\n",options->cloud_points_filename);
    else
        printf("cloud_points_filename = NULL\n");
    printf("********************************************************\n");
}