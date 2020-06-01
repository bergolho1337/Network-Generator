#include "network/tree.h"
#include "options/user_options.h"

using namespace std;

int main (int argc, char *argv[])
{
	if (argc-1 != 1)
	{
		usage(argv[0]);
		exit(EXIT_FAILURE);
	}

	User_Options *the_options = new User_Options(argc,argv);
	Tree *the_tree = new Tree(the_options);

	the_tree->generate();

	delete the_tree;
	delete the_options;

	return 0;
}
