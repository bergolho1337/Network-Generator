#include "../include/phasefield.h"

PhaseField_Config::PhaseField_Config () { }

void PhaseField_Config::print ()
{
	printf("[PhaseField] root_point = (%.10lf,%.10lf,%.10lf)\n",this->root_point[0],this->root_point[1],this->root_point[2]);
	printf("[PhaseField] miocardium_filename = %s\n",this->miocardium_filename.c_str());
	printf("[PhaseField] terminals_filename = %s\n",this->terminals_filename.c_str());
}

PhaseField_Generator::PhaseField_Generator (PhaseField_Config *config)
{
	printf("[PhaseField] Building PhaseField network ...\n");
	printf("[PhaseField] Need to be implemented yet ...\n");	

	// Here comes the PhaseField algorithm to build the network
	// ...
}

void PhaseField_Generator::write_network_to_VTK ()
{
	// Write the network into a VTK file ...
}
