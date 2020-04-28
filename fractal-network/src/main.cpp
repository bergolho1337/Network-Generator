#include <cstdio>
#include <cstdlib>

#include "options/user_options.h"
#include "fractal/fractal.h"

using namespace std;

void configure_fractal_tree_from_options (Fractal_Tree *the_tree, User_Options *the_options)
{
	the_tree->max_iterations = the_options->max_iterations;
	the_tree->initial_length = the_options->initial_length;
	the_tree->initial_angle = the_options->initial_angle;
	the_tree->initial_diameter = the_options->initial_diameter;
	the_tree->angle_decrease_ratio = the_options->angle_decrease_ratio;
	the_tree->length_decrease_ratio = the_options->length_decrease_ratio;
	the_tree->diameter_decrease_ratio = the_options->diameter_decrease_ratio;
	for (uint32_t i = 0; i < 3; i++)
		the_tree->root_pos[i] = the_options->root_pos[i];
}

int main (int argc, char *argv[])
{
	if (argc-1 != 1)
	{
		usage(argv[0]);
		exit(EXIT_FAILURE);
	}

	User_Options *the_options = new User_Options(argc,argv);
	Fractal_Tree *the_tree = new Fractal_Tree();

	configure_fractal_tree_from_options(the_tree,the_options);

	the_tree->grow_network();
	the_tree->write();
	//the_tree->print();

	delete the_tree;
	delete the_options;

	return 0;
}
