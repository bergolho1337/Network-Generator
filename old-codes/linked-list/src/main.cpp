// Author: Lucas Berg

#include <iostream>
#include <cstdio>
#include <cstdlib>

#include "linked-list/linked-list.h"

using namespace std;

int main (int argc, char *argv[])
{   
    // Test 1
    struct linked_list *list1 = new_linked_list();
    insert_node(list1,5.0);
    insert_node(list1,3.0);
    insert_node(list1,10.0);
    insert_node(list1,2.0);
    print_list(list1);

    delete_node(list1,10.0);
    print_list(list1);

    free_linked_list(list1);

    // Test 2
    struct linked_list *list2 = new_linked_list();
    for (int i = 0; i < 100; i++)
        insert_node(list2,i);
    print_list(list2);
    free_linked_list(list2);

    // Test 3
    struct linked_list *list3 = new_linked_list();
    insert_node(list3,1.0);
    insert_node(list3,2.0);
    insert_node(list3,3.0);
    struct node *tmp1 = search_node(list3,2.0);
    struct node *tmp2 = search_node(list3,100.0);
    free_linked_list(list3);

    // Test 4
    struct linked_list *list4 = new_linked_list();
    insert_node(list4,10.0);
    insert_node(list4,2.0);
    insert_node(list4,7.0);
    insert_node(list4,1.0);
    insert_node(list4,5.0);
    print_list(list4);

    delete_node(list4,1.0);
    order_list(list4);
    print_list(list4);

    free_linked_list(list4);


    return 0;
}