#ifndef GENERATOR_H
#define GENERATOR_H

#include <iostream>
#include <fstream>
#include <string>
#include <cstdio>
#include <cstdlib>
#include <cmath>

#include "config.h"

class Network_Generator
{
private:
	// Network generator algorithm configuration
	std::string method_name;

	// Structures for the generation algorithm
	// TODO: Make this an abstract class
	Lsystem_Generator *lsystem_generator;
	CO_Generator *co_generator;
	PhaseField_Generator *phase_generator;

public:
	Network_Generator (User_Options *options);
	~Network_Generator ();
};

#endif
