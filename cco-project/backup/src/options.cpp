#include "../include/options.h"

User_Options* read_user_input (int argc, char *argv[])
{
    User_Options *options = new User_Options();

    options->Q_perf = atof(argv[1]);
    options->p_perf = atof(argv[2]);
    options->r_perf = atof(argv[3]);
    options->N_term = atoi(argv[4]);

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