#include "../include/generator.h"

Network_Generator::Network_Generator (User_Options *options)
{
	this->method_name = options->method_name;
	if (this->method_name == "Lsystem")
	{
		lsystem_generator = new Lsystem_Generator(options->lsystem_config);
		lsystem_generator->write_network_to_VTK();

		co_generator = NULL;
		phase_generator = NULL;
	}
	else if (this->method_name == "CO")
	{
		co_generator = new CO_Generator(options->co_config);
		co_generator->write_network_to_VTK();

		lsystem_generator = NULL;
		phase_generator = NULL;
	}	
	else if (this->method_name == "PhaseField")
	{
		phase_generator = new PhaseField_Generator(options->phase_config);
		phase_generator->write_network_to_VTK();

		lsystem_generator = NULL;
		co_generator = NULL;
	}
	else
	{
		printf("[-] ERROR! Invalid method name\n");
		exit(EXIT_FAILURE);
	}
}

Network_Generator::~Network_Generator ()
{
	if (lsystem_generator)
		delete lsystem_generator;

	//if (co_generator)
	//	delete co_generator;

	//if (phase_generator)
	//	delete phase_generator;
}

