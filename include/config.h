#ifndef CONFIG_H
#define CONFIG_H

#include <iostream>
#include <fstream>
#include <string>
#include <cstdio>
#include <cstdlib>
#include <cmath>

#include "lsystem.h"
#include "co.h"
#include "phasefield.h"

class User_Options
{
public:
	// Network generator algorithm configuration
	std::string method_name;

	// Reference to the method object
	// TODO: Change this to an abstract class ...
	Lsystem_Config *lsystem_config;
	CO_Config *co_config;
	PhaseField_Config *phase_config;

	// Private functions to the class
	void Read_Input_File (const char filename[]);
public:
	User_Options (int argc, char *argv[]);
	void print ();
	
};

#endif
