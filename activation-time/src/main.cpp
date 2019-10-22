// Author: Lucas Berg

#include <iostream>
#include <vector>
#include <algorithm>

#include <cstdio>
#include <cstdlib>
#include <cmath>

using namespace std;

// ====================================================================================================================
// CONSTANTS AND MACROS
#define PRINT_LINE "--------------------------------------------------------------------------------------------"
#define UM_TO_CM 0.0001
#define MS_TO_S 0.001
#define S_TO_MS 1000.0
#define CM_S_TO_M_S 0.01

const double C = 1080.8;
const double CM = 1.2;                      // uF/cm^2
const double SIGMA = 0.004;                 // mS/cm
const double RM = 6750.0;                   // ohm.cm^2
const double RC = 150.0;                    // ohm.cm

const double SIGMA_M = 0.15;                // cm
const double TAU_M = 8.4;                   // ms

const double MIN_RADIUS = 10.0;             // um
const double MAX_RADIUS = 3000.0;             // um
const double CABLE_LENGTH = 1.0;            // cm
const unsigned int NPOINTS = 500;
// ====================================================================================================================

// Output will be given in (cm)
double calc_lambda_m (const double r)
{
    double d = 2.0 * r;

    // Equation (4.17) from Keener's book
    double num = RM * d * UM_TO_CM;
    double den = 4.0 * RC;

    return sqrt(num / den);
}

// Output will be given in (s)
double calc_tau_m ()
{
    // Equation (4.16) from Keener's book
    return RM * CM * MS_TO_S;
}

// Output will be given in (cm/s)
double calc_velocity (const double r)
{
    return C * calc_lambda_m(r) / calc_tau_m();
}

// Ouptut will be given in (ms)
double calc_activation_time (const double r)
{
    return CABLE_LENGTH / calc_velocity(r) * S_TO_MS;
}

void write_datafile ()
{
    FILE *file_lambda = fopen("output/lambda.dat","w+");
    FILE *file_velocity = fopen("output/velocity.dat","w+");
    FILE *file_at = fopen("output/activation_time.dat","w+");

    double delta = (MAX_RADIUS - MIN_RADIUS) / (double)NPOINTS;

    for (unsigned int i = 0; i < NPOINTS; i++)
    {
        double radius = MIN_RADIUS + delta*i;

        double lambda_m = calc_lambda_m(radius);
        double velocity = calc_velocity(radius);
        double at = calc_activation_time(radius);

        fprintf(file_lambda,"%g %g\n",radius,lambda_m);
        fprintf(file_velocity,"%g %g\n",radius,velocity);
        fprintf(file_at,"%g %g\n",radius,at);
    }

    fclose(file_lambda);
    fclose(file_velocity);
    fclose(file_at);
}

void write_function ()
{
    FILE *file_new_at = fopen("output/new_cost_function.dat","w+");

    double dx = 0.1;
    double base = 2.0;

    for (int i = 0; i < 100; i++)
    {
        double x = i*dx;

        double y = pow(base,-1/x);

        fprintf(file_new_at,"%g %g\n",x,y);
    }

    fclose(file_new_at);
}

int main (int argc, char *argv[])
{
    if (argc-1 != 0)
    {
        printf("%s\n",PRINT_LINE);
        printf("Usage:> %s\n",argv[0]);
        printf("%s\n",PRINT_LINE);
        exit(EXIT_FAILURE);
    }

    write_datafile();
    write_function();

    return 0;
}
