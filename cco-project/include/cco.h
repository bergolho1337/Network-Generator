#ifndef CCO_H
#define CCO_H

#include <ctime>

#include "graph.h"
#include "options.h"

// Constant and macros
const double GAMMA = 2.55;       // Bifurcation expoent
const double ETA = 3.6e-03;      // Viscosity of blood

bool is_inside_circle (const double x, const double y, const double r);
void make_root (Graph *the_network, const double Q_perf, const double Nterm);
void grow_cco_tree (Graph *the_network, User_Options *options);

#endif