#include "user_options.h"

struct user_options* new_user_options (int argc, char *argv[])
{
    struct user_options *result = (struct user_options*)malloc(sizeof(struct user_options));
    result->Q_perf = atof(argv[1]);
    result->p_perf = atof(argv[2]);
    result->p_term = atof(argv[3]);
    result->r_perf = atof(argv[4]);
    result->N_term = atoi(argv[5]);
    if (argc-1 == 6)
        strcpy(result->cloud_filename,argv[6]);
    else
        strcpy(result->cloud_filename,"");
    
    return result;
}

void free_user_options (struct user_options *options)
{
    free(options);
}