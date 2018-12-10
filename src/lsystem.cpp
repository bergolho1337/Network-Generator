#include "../include/lsystem.h"

// Build a Lsystem Purkinje network over a endocardium mesh
Lsystem_Generator::Lsystem_Generator (Lsystem_Config *config)
{
	printf("[Lsystem] Building Lsystem network ...\n");

	this->the_purkinje_network = new Graph();
	this->the_miocardium = new Miocardium();

	this->the_miocardium->read_cloud_points(config->miocardium_filename);
	this->the_miocardium->read_terminal_points(config->terminals_filename);

	this->the_miocardium->set_limits();	

	this->initialize_root_point();
	this->initialize_random_array(config->branch_length);
	
	// TODO: Implement the rest of the functions ...

}



Lsystem_Config::Lsystem_Config ()
{

}

// TODO: Put this option in the input configuration file
void Lsystem_Generator::initialize_root_point ()
{
	// Benjamin's root position : (158.5,240.28,54.6504)
	// Scaling to the reduced mesh ...
	this->root_point[0] = 39.625;
	this->root_point[1] = 60.07;
	this->root_point[2] = 13.6626;
}

void Lsystem_Generator::initialize_random_array (const double l_bra)
{
	srand(time(NULL));
	default_random_engine generator;
	normal_distribution<double> distribution(0.0,l_bra*SIGMA_RANDOM_DISTRIBUTION);			
	for (int i = 0; i < MAX_SIZE_RAND_ARRAY; i++)
		this->rand_numbers[i] = distribution(generator);
}


void Lsystem_Config::print ()
{
	printf("[Lsystem] max_grow_iterations = %u\n",this->max_grow_iterations);
	printf("[Lsystem] branch_length = %.10lf\n",this->branch_length);
	printf("[Lsystem] w_1 = %.10lf\n",this->w_1);
	printf("[Lsystem] tolerance_tree_collision = %.10lf\n",this->tolerance_tree_collision);
	printf("[Lsystem] tolerance_miocardium_collision = %.10lf\n",this->tolerance_miocardium_collision);
	printf("[Lsystem] tolerance_terminal_collision = %.10lf\n",this->tolerance_terminal_collision);
	printf("[Lsystem] miocardium_filename = %s\n",this->miocardium_filename.c_str());
	printf("[Lsystem] terminals_filename = %s\n",this->terminals_filename.c_str());
	
}

