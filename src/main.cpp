// --------------------------------------------------------------------------------------------------------------------
// This program builds several types of Purkinje network using different techniques:
// 	L-System
//	CCO
//	Phase Field
// --------------------------------------------------------------------------------------------------------------------

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <omp.h>

#include "../include/generator.h"
#include "../include/timer.h"
#include "../include/config.h"
#include "../include/utils.h"
//#include "../include/plot.h"
//#include "../include/stimuli.h"
//#include "../include/monodomain.h"

using namespace std;

int main (int argc, char *argv[])
{
	if (argc-1 != 1)
	{
		Usage(argv[0]);
		exit(EXIT_FAILURE);
	}

	// Initialize all the structures 
	User_Options *options = new User_Options(argc,argv);

	Network_Generator *generator = new Network_Generator(options);

	/*	
	struct monodomain_solver *solver;
	solver = new_monodomain_solver();

	struct plot_config *plotter;
	plotter = new_plot_config();

	struct stim_config *stim;
	stim = new_stim_config();

	// Reading input configuration file
	read_input_file(options,argv[1]);

	// Parse the user input into the structures
	configure_solver_from_options(solver,options);
	configure_plot_from_options(plotter,solver,options);
	configure_stimulus_from_options(stim,options);

	// OpenMP configuration
	omp_set_dynamic(0);
    omp_set_num_threads(4);

	print_configuration_parameters(solver->dx,solver->dt,solver->tmax,solver->lmax,\
					solver->Ncell,solver->Niter,Nodes,\
					plotter->plot_cell_ids,plotter->print_rate,plotter->sst_rate);

	// Call the solver function for the monodomain equation
	solve_monodomain(solver,stim,plotter);
	
	// Free memory
	free_user_options(options);
	free_monodomain_solver(solver);
	free_plot_config(plotter);
	free_stim_config(stim);
	*/

	return 0;  
}
