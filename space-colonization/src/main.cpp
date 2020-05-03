#include "network/tree.h"

using namespace std;

int main (int argc, char *argv[])
{
	Tree *the_tree = new Tree();

	for (uint32_t i = 0; i < 200; i++)
	{
		the_tree->grow_network();
		the_tree->write_branches(i);
		the_tree->write_leaves(i);
	}

	delete the_tree;

	return 0;
}
