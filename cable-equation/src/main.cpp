// Author: Lucas Berg
// Most of the information to adjust the celular parameter 'c' came from two sources:
//  KEENER, James P.; SNEYD, James. Mathematical physiology. New York: Springer, 1998.
//  FOZZARD, Harry A. Membrane capacity of the cardiac Purkinje fibre. The Journal of physiology, v. 182, n. 2, p. 255-267, 1966.

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
const double MAX_RADIUS = 1500.0;             // um
const double CABLE_LENGTH = 5.0;            // cm
const unsigned int NPOINTS = 500;
// ====================================================================================================================

// Output will be given in (cm)
double calc_lambda_m (const double Rm, const double Rc, const double d)
{

    // Equation (4.17) from Keener's book
    double num = 1000.0 * Rm * d * UM_TO_CM;         // Convert the diameter to {cm}
    double den = 4.0 * Rc;

    return sqrt(num / den);
}

// Output will be given in (ms)
double calc_tau_m (const double Rm, const double Cm)
{
    // Equation (4.16) from Keener's book
    return Rm * Cm;
}

// Output will be given in (cm/s)
double calc_velocity (const double lambda, const double tau, const double c)
{
    return c * lambda / tau;
}

/*
// Ouptut will be given in (ms)
double calc_activation_time (const double r)
{
    return CABLE_LENGTH / calc_velocity(r) * S_TO_MS;
}
*/

void write_datafile (const double Rm, const double Rc, const double Cm, const double c)
{
  FILE *file_s = fopen("output/propagation_velocity.txt","w+");
  FILE *file_at = fopen("output/activation_time.txt","w+");

  double delta = (MAX_RADIUS - MIN_RADIUS) / (double)NPOINTS;

  for (uint32_t i = 0; i < NPOINTS; i++)
  {
    double radius = MIN_RADIUS + delta*i;
    double diameter = 2.0*radius;

    double lambda_m = calc_lambda_m(Rm,Rc,diameter);
    double tau_m = calc_tau_m(Rm,Cm);
    double s = calc_velocity(lambda_m,tau_m,c)*10.0;    // mm/ms
    double at = CABLE_LENGTH / s;                       // ms

    fprintf(file_s,"%g %g\n",diameter,s);
    fprintf(file_at,"%g %g\n",diameter,at);
  }
  fclose(file_s);
  fclose(file_at);
}

int main (int argc, char *argv[])
{
    if (argc-1 < 4)
    {
        printf("%s\n",PRINT_LINE);
        printf("Usage:> %s <Rm> <Rc> <Cm> <d> [s] [lambda_m] [tau_m]\n",argv[0]);
        printf("%s\n",PRINT_LINE);
        printf("<Rm> = Membrane resistance {k.ohm.cm^2}\n");
        printf("<Rc> = Citoplasm resistance {ohm.cm}\n");
        printf("<Cm> = Membrane capacitance {uF/cm^2}\n");
        printf("<d> = Diameter {um}\n");
        printf("[s] = Experimental propagation velocity {mm/ms}\n");
        printf("%s\n",PRINT_LINE);
        exit(EXIT_FAILURE);
    }

    // Purkinje muscle cell (mean) --> human
    // Rm = 1714 ohm.cm^2
    // Rc = 116 ohm.cm
    // Cm = 12.8 uF/cm^2
    // s = 2.667 m/s
    // d = 78 um
    // c = 28.21

    // Purkinje muscle cell (mean) --> Fenton experiment
    // Rm = 1714 ohm.cm^2
    // Rc = 116 ohm.cm
    // Cm = 12.8 uF/cm^2
    // s = [2,4] m/s
    // d = [2000,3000] um
    // c = 8.51

    //const double c = 8.51;        // Fenton experiment (Canine Purkinje cell)
    const double c = 28.21;       // Human Purkinje cell
    //const double c = 3.87298;     // Good aproximation for the propagation velocity on the majority of the excitable cells

    // Parse the input
    double Rm = atof(argv[1]);
    double Rc = atof(argv[2]);
    double Cm = atof(argv[3]);
    double d = atof(argv[4]);

    double lambda_m = calc_lambda_m(Rm,Rc,d);
    double tau_m = calc_tau_m(Rm,Cm);

    printf("Lambda_m = %.2lf cm\n",lambda_m);
    printf("Tau_m = %.2lf ms\n",tau_m);
    printf("Calculated velocity = %g mm/ms\n",calc_velocity(lambda_m,tau_m,c)*10.0);
    if (argc-1 == 7)
    {
      double s = atof(argv[5]);
      double lambda_m_exp = atof(argv[6]);
      double tau_m_exp = atof(argv[7]);
      double c_exp = (s*0.1*tau_m_exp)/(lambda_m_exp);
      printf("Experiment lambda_m = %g cm\n",lambda_m_exp);
      printf("Experiment tau_m = %g ms\n",tau_m_exp);
      printf("Experiment velocity = %g mm/ms\n",s);
      printf("Experiment c value = %g\n",c_exp);

    }

    write_datafile(Rm,Rc,Cm,c);

    return 0;
}
