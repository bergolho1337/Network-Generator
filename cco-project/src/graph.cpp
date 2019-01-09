#include "../include/graph.h"

Graph* initialize_graph ()
{
    Graph *the_network = new Graph();

    the_network->list_nodes = NULL;
    the_network->last_node = NULL;
    the_network->num_nodes = 0;
    the_network->num_segments = 0;

    return the_network;
}

void insert_node_graph (Graph *g, const double x, const double y)
{
    Node *new_node = new Node();
    new_node->index = g->num_nodes;
    new_node->x = x;
    new_node->y = y;
    new_node->next = NULL;
    new_node->parent = NULL;
    new_node->left_offspring = NULL;
    new_node->right_offpring = NULL;
    new_node->segment_list = NULL;
    new_node->num_segments = 0;

    if (g->last_node == NULL)
    {
        g->list_nodes = new_node;
        g->last_node = new_node;
        g->num_nodes++;
    }
    else
    {
        g->last_node->next = new_node;
        g->last_node = new_node;
        g->num_nodes++;
    }
}

Node* search_node (Graph *g, const int id)
{
    Node *ptr = g->list_nodes;
    while (ptr != NULL)
    {
        if (ptr->index == id)
            return ptr;
        ptr = ptr->next;
    }
    printf("[-] ERROR! Node %d not found!\n",id);
    return NULL;
}

double calc_segment_length (const Node *ptr1, const Node *ptr2)
{
    return sqrt(pow(ptr1->x-ptr2->x,2) + pow(ptr1->y-ptr2->y,2));
}

void insert_edge_graph (Graph *g, const int source, const int destination)
{
    Node *ptr1 = search_node(g,source);
    Node *ptr2 = search_node(g,destination);

    Edge *new_segment = new Edge();
    new_segment->index = destination;
    new_segment->length = calc_segment_length(ptr1,ptr2);
    //new_segment->radius = poisseule(length);
    //new_segment->resistance = ???
    //new_segment->pressure = ???
    new_segment->next = NULL;

    Edge *ptrl = ptr1->segment_list;
    if (ptrl == NULL)
    {
        ptr1->segment_list = new_segment;
    }
    else
    {
        Edge *ptrl = ptr1->segment_list;
        while (ptrl->next != NULL)
        {
            ptrl = ptrl->next;
        }
        ptrl->next = new_segment;
    }
    ptr1->num_segments++;
    g->num_segments++;

    // TODO
    // Adjust pointer to parent and offsprings ..
}

void print_graph (Graph *g)
{
    Node *ptr = g->list_nodes;
    while (ptr != NULL)
    {
        printf("|| %d (%lf %lf) ||", ptr->index,ptr->x,ptr->y);
        Edge *ptrl = ptr->segment_list;
        while (ptrl != NULL)
        {
            printf(" --> || %d %lf ||",ptrl->index,ptrl->length);
            ptrl = ptrl->next;
        }
        printf("\n");
        ptr = ptr->next;
    }
}

void write_graph_to_VTK (Graph *g)
{
    Node *ptr;
    Edge *ptrl;

    FILE *file = fopen("cco_tree.vtk","w+");

    fprintf(file,"# vtk DataFile Version 3.0\n");
    fprintf(file,"CCO\n");
    fprintf(file,"ASCII\n");
    fprintf(file,"DATASET POLYDATA\n");
    fprintf(file,"POINTS %d float\n",g->num_nodes);
    ptr = g->list_nodes;
    while (ptr != NULL)
    {
        fprintf(file,"%g %g 0\n",ptr->x,ptr->y);
        ptr = ptr->next;
    }
    fprintf(file,"LINES %d %d\n",g->num_segments,g->num_segments*3);
    ptr = g->list_nodes;
    while (ptr != NULL)
    {
        ptrl = ptr->segment_list;
        while (ptrl != NULL)
        {
            fprintf(file,"2 %d %d\n",ptr->index,ptrl->index);
            ptrl = ptrl->next;
        }
        ptr = ptr->next;
    }

    fclose(file);
}