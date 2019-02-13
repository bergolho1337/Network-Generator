// Author: Lucas Berg

#include <iostream>
#include <cstdio>
#include <cstdlib>

//#include "point-list/point-list.h"
//#include "segment-list/segment-list.h"
#include "cco/cco.h"
#include "options/user_options.h"

using namespace std;

int main (int argc, char *argv[])
{   
    if (argc-1 != 4)
    {
        printf("Usage:> %s <Qperf> <pperf> <rperf> <Nterm>\n",argv[0]);
        exit(EXIT_FAILURE);   
    }

    struct user_options *options = new_user_options(argc,argv);
    struct cco_network *the_network = new_cco_network(options);

    grow_tree(the_network);
    write_to_vtk(the_network);

    /*
    // Test Point list
    struct point_list *list1 = new_point_list();
    double pos1[3] = {0,0,0};
    double pos2[3] = {1,1,0};
    double pos3[3] = {2,2,0};
    double pos4[3] = {3,3,0};
    
    insert_node(list1,pos1);
    insert_node(list1,pos2);
    insert_node(list1,pos3);
    insert_node(list1,pos4);
    
    print_list(list1);

    free_point_list(list1);

    // Test Segment list
    struct segment_list *list2 = new_segment_list();
    struct segment *seg1 = new_segment(0,1,NULL,NULL,NULL,1,1);
    struct segment *seg2 = new_segment(1,2,NULL,NULL,NULL,1,1);
    struct segment *seg3 = new_segment(1,3,NULL,NULL,NULL,1,1);
    insert_node(list2,seg1);
    insert_node(list2,seg2);
    insert_node(list2,seg3);

    print_list(list2);

    free_segment_list(list2);
    */

    return 0;
}