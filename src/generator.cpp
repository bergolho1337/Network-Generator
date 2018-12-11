#include "../include/generator.h"

Network_Generator::Network_Generator (User_Options *options)
{
	this->method_name = options->method_name;
	if (this->method_name == "Lsystem")
	{
		lsystem_generator = new Lsystem_Generator(options->lsystem_config);
		lsystem_generator->write_network_to_VTK();
	}
	/*
	else if (this->method_name == "CO")
		co_generator = new CO_Generator(options->co_config);
	if (this->method_name == "PhaseField")
		phase_generator = new PhaseField_Generator(options->phase_config);
	*/

}
