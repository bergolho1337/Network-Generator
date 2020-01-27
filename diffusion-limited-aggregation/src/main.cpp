// Author: Lucas Berg

#include <iostream>
#include <cstdio>
#include <cstdlib>

#include "tree/tree.h"

using namespace std;

int main (int argc, char *argv[])
{   
    Tree *the_tree = new Tree();

    the_tree->grow();

    return 0;
}