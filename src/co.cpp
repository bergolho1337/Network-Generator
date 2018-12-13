#include "../include/co.h"

CO_Config::CO_Config () { }

void CO_Config::print ()
{
	printf("[CO] root_point = (%.10lf,%.10lf,%.10lf)\n",this->root_point[0],this->root_point[1],this->root_point[2]);
	printf("[CO] miocardium_filename = %s\n",this->miocardium_filename.c_str());
	printf("[CO] terminals_filename = %s\n",this->terminals_filename.c_str());
}

CO_Generator::CO_Generator (CO_Config *config)
{
	printf("[CO] Building CO network ...\n");
	printf("[CO] Need to be implemented yet ...\n");	

	// Here comes the CO algorithm to build the network
	// ...
}

void CO_Generator::write_network_to_VTK ()
{
	// Write the network into a VTK file ...
}
