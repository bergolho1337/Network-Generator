// Author: Lucas Berg

#include <iostream>
#include <cstdio>
#include <cstdlib>

#include "tree/tree.h"
#include "options/user_options.h"

using namespace std;

int main (int argc, char *argv[])
{   
    if (argc-1 != 1)
    {
        usage(argv[0]);
        exit(EXIT_FAILURE);
    }

    struct user_options *the_options = new_user_options(argc,argv);
    //print_user_options(the_options);

    struct dla_tree *the_tree = new_dla_tree();
    //print_dla_tree(the_tree);

    grow_tree(the_tree, the_options);
    //print_dla_tree(the_tree);

    free_dla_tree(the_tree);
    free_user_options(the_options);

    return 0;
}