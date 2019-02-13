#ifndef OPTIONS_H
#define OPTIONS_H

#include <cstdio>
#include <cstdlib>
#include <cstring>

struct User_Options
{
    double x0, y0;              // Proximal point of the root
    double Q_perf;              // Perfusion flow of the root
    double p_perf;              // Perfusion pressure of the root
    double r_perf;              // Perfusion radius
    int N_term;                 // Number of terminals

}typedef User_Options;

User_Options* read_user_input (int argc, char *argv[]);
void print_user_input (User_Options *options);

#endif