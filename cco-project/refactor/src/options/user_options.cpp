#include "user_options.h"

struct user_options* new_user_options (int argc, char *argv[])
{
    struct user_options *result = (struct user_options*)malloc(sizeof(struct user_options));
    result->Q_perf = atof(argv[1]);
    result->p_perf = atof(argv[2]);
    result->r_perf = atof(argv[3]);
    result->N_term = atoi(argv[4]);
    return result;
}