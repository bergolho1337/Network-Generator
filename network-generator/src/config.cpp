#include "../include/config.h"

User_Options::User_Options (int argc, char *argv[])
{
	printf("[User_Options] Reading user input ...\n");
	Read_Input_File(argv[1]);
}

User_Options::~User_Options ()
{
	
}

void User_Options::Read_Input_File (const char filename[])
{
	std::string str, str2, str3;
	std::ifstream file(filename);
		
	// Reading main section
	while (getline(file,str) && str != "[main]");
	file >> str >> str2 >> method_name;
	
	// Switch the reading format using the method_name
	if (method_name == "Lsystem")
	{
		printf("[User_Options] Setting Lsystem configuration ...\n");
		lsystem_config = new Lsystem_Config();
		file >> str >> str2 >> lsystem_config->root_point[0]; 
		file >> str >> str2 >> lsystem_config->root_point[1];
		file >> str >> str2 >> lsystem_config->root_point[2];
		file >> str >> str2 >> lsystem_config->max_grow_iterations;
		file >> str >> str2 >> lsystem_config->branch_length;
		file >> str >> str2 >> lsystem_config->w_1;
		file >> str >> str2 >> lsystem_config->sobel_cube_size;
		file >> str >> str2 >> lsystem_config->tolerance_tree_collision;
		file >> str >> str2 >> lsystem_config->tolerance_miocardium_collision;
		file >> str >> str2 >> lsystem_config->tolerance_terminal_collision;
		file >> str >> str2 >> lsystem_config->miocardium_filename;
		file >> str >> str2 >> lsystem_config->terminals_filename;

		// DEBUG
		//lsystem_config->print();

	}
	else if (method_name == "CO")
	{
		printf("[User_Options] Setting CO configuration ...\n");
		co_config = new CO_Config();
		file >> str >> str2 >> co_config->root_point[0]; 
		file >> str >> str2 >> co_config->root_point[1];
		file >> str >> str2 >> co_config->root_point[2];
		file >> str >> str2 >> co_config->miocardium_filename;
		file >> str >> str2 >> co_config->terminals_filename;

		// DEBUG
		//co_config->print();
	}
	else if (method_name == "PhaseField")
	{
		printf("[User_Options] Setting PhaseField configuration ...\n");
		phase_config = new PhaseField_Config();
		file >> str >> str2 >> phase_config->root_point[0]; 
		file >> str >> str2 >> phase_config->root_point[1];
		file >> str >> str2 >> phase_config->root_point[2];
		file >> str >> str2 >> phase_config->miocardium_filename;
		file >> str >> str2 >> phase_config->terminals_filename;

		// DEBUG
		//phase_config->print();
	}
	else
	{
		printf("[User_Options] ERROR! Invalid method name! '%s' does not name a method!\n",method_name.c_str());
		exit(EXIT_FAILURE);
	}
	
	file.close();
}

void User_Options::print ()
{
	printf("Method Name = %s\n",method_name.c_str());
	lsystem_config->print();
}

