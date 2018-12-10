#ifndef UTILS_H
#define UTILS_H

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#include <cstdio>
#include <cstdlib>
#include <cmath>

#define PRINT_LINE "===================================================================================================="
/*
void write_VTK_to_file (const double *sv, const double dx,\
                        const int Ncell, const int Nodes,\
			const int iter);
void write_plot_data (std::ofstream files[], const double t, const double *sv,\
			const int Ncell, const int Nodes,\
			const int ids[]);
void write_steady_state_to_file (const double *sv, const int Ncell, const int Nodes);

void print_stimulus (const double *stim_current, const int Ncell, const double dx);
void print_state_vector (const double *sv, const int Ncell, const int Nodes);
void print_progress (int iter, int max_iter);
void print_configuration_parameters(const double dx, const double dt, const double tmax, const double lmax,\
					const int Ncell, const int Niter, const int Nodes,\
					const int plot_cell_ids[], const int print_rate, const int sst_rate);

void configure_plot_cells (std::ofstream files[], int plot_cell_ids[], const double lmax, const double dx);
void close_plot_files (std::ofstream files[]);
*/
void Usage (const char pname[]);


#endif
