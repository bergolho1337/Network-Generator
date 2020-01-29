// Author: Lucas Berg

#include <iostream>
#include <cstdio>
#include <cstdlib>

#include "tree/tree.h"
#include "options/user_options.h"

using namespace std;

/*
// TODO: Move this to another folder ...
void grow_tree (DLA_Tree *the_tree, User_Options *the_options)
{
    // Make the root
    Walker *root = new Walker(WIDTH/2.0,HEIGHT/2.0,0.0);
    the_tree->the_points.push_back(root);

    // Add the Walkers
    std::vector<Walker*> the_others;
    for (uint32_t i = 0; i < MAX_NUMBER_OF_WALKERS; i++)
    {
        Walker *walker = new Walker();
        //this->the_tree.push_back(walker);
        the_others.push_back(walker);
    }

    // Main iteration loop 
    for (uint32_t iter = 0; iter < MAX_NUMBER_OF_ITERATIONS; iter++)
    {
        print_progress_bar(iter,MAX_NUMBER_OF_ITERATIONS);

        // Move each Walker using the Random Walk
        for (uint32_t i = the_others.size()-1; i > 0; i--)
        {
            //printf("%u\n",i);
            the_others[i]->walk();

            uint32_t stuck_index = the_others[i]->is_stuck(the_tree->the_points);
            if (stuck_index != the_tree->the_points.size())
            {
                uint32_t new_index = the_tree->the_points.size();
                Segment *segment = new Segment(stuck_index,new_index);
                the_tree->the_segments.push_back(segment);

                the_tree->the_points.push_back(the_others[i]);
                the_others.erase(the_others.begin()+i);

                the_tree->write_to_vtk(the_tree->the_points.size());
            }
        }

        // Add more walkers as the tree grows ...
        while (the_others.size() < MAX_NUMBER_OF_WALKERS)
        {
            Walker *walker = new Walker();
            the_others.push_back(walker);
        }
    }
    printf("\n");
}
*/

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
    grow_tree(the_tree, the_options);

    free_dla_tree(the_tree);
    free_user_options(the_options);

    return 0;
}