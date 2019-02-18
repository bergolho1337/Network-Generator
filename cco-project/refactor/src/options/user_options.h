//
// Created by bergolho on 12/02/19.
//

#ifndef USER_OPTIONS_H
#define USER_OPTIONS_H

#include <cstdio>
#include <cstdlib>
#include <cstdbool>
#include <cstdint>
#include "../point-list/point-list.h"

struct user_options
{
    int N_term;
    double Q_perf;
    double p_perf;
    double p_term;
    double r_perf;
};

struct user_options* new_user_options (int argc, char *argv[]);
void free_user_options (struct user_options *options);

#endif