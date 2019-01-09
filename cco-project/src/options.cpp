#include "../include/options.h"

User_Options* read_user_input (int argc, char *argv[])
{
    User_Options *options = new User_Options();

    options->x0 = atof(argv[1]); options->y0 = atof(argv[2]);
    options->Q_perf = atof(argv[3]);
    options->p_perf = atof(argv[4]);
    options->r_perf = atof(argv[5]);
    options->N_term = atoi(argv[6]);

    return options;
}

void print_user_input (User_Options *options)
{
    printf("Root proximal point:\n");
    printf("\t(%.10lf %.10lf)\n",options->x0,options->y0);
    printf("Qperf = %.10lf\n",options->Q_perf);
    printf("pperf = %.10lf\n",options->p_perf);
    printf("rperf = %.10lf\n",options->r_perf);
    printf("Nterm = %d\n",options->N_term);
}