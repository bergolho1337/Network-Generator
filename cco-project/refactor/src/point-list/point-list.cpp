#include "point-list.h"

struct point_list* new_point_list ()
{
    struct point_list *l = (struct point_list*)malloc(sizeof(struct point_list));
    l->num_nodes = 0;
    l->list_nodes = NULL;
    return l;
}

void free_point_list (struct point_list *l)
{
    uint32_t cont = 0;
    while (!is_empty(l))
    {
        delete_node(l,cont);
        cont++;
    }
    //print_list(l);
    
    free(l);
    l = NULL;
}

struct point_node* insert_point (struct point_list *l, const double pos[])
{
    struct point_node *node;
    // List is empty
    if (!l->list_nodes)
    {
        node = new_node(l->num_nodes,pos);
        l->list_nodes = node;
    }
    // List has some values
    else
    {
        // Pass through all the list
        struct point_node *tmp = l->list_nodes;
        while (tmp->next != NULL)
        {
            tmp = tmp->next;
        }

        node = new_node(l->num_nodes,pos);
        tmp->next = node;
    }
    l->num_nodes++;
    return node;
}

void delete_node (struct point_list *l, const uint32_t index)
{
    struct point_node *aux1, *aux2;
    aux1 = l->list_nodes;
    aux2 = NULL;

    if (l->num_nodes == 0)
    {
        fprintf(stderr,"[-] ERROR! The list is empty!\n");
        exit(EXIT_FAILURE);
    }

    while (aux1 != NULL && aux1->id != index)
    {
        aux2 = aux1;
        aux1 = aux1->next;
    }

    // Case 1: First element
    if (aux1 == l->list_nodes)
    {
        //printf("Removing %lf --> Case 1\n",key);
        l->list_nodes = aux1->next;
        aux1->next = NULL;
    }
    // Case 2: Last element
    else if (aux1->next == NULL)
    {
        //printf("Removing %lf --> Case 2\n",key);
        aux2->next = NULL;
    }
    // Case 3: Middle element
    else
    {
        //printf("Removing %lf --> Case 3\n",key);
        aux2->next = aux1->next;
        aux1->next = NULL;
    }
    free_point(aux1->value);
    free(aux1);
    l->num_nodes--;
}

struct point_node* search_node (struct point_list *l, const uint32_t index)
{
    struct point_node *tmp = l->list_nodes;
    while (tmp != NULL)
    {
        if (tmp->id == index)
        {
            printf("[+] Found the node || %u %lf ||\n",tmp->id,tmp->value);
            return tmp;
        }
            
        tmp = tmp->next;
    }
    fprintf(stderr,"[-] Index %u not found!\n",index);
    return NULL;
}

bool is_empty (struct point_list *l)
{
    return (l->list_nodes == NULL) ? true : false;
}

void print_list (struct point_list *l)
{
    if (l->num_nodes == 0)
    {
        printf("[!] The list is empty!\n");
        return;
    }

    struct point_node *tmp = l->list_nodes;
    printf("Number of node in list = %u\n",l->num_nodes);
    while (tmp != NULL)
    {
        printf("Node %d -- Point(%g,%g,%g)\n",tmp->id,\
                                        tmp->value->x,\
                                        tmp->value->y,\
                                        tmp->value->z);
        tmp = tmp->next;
    }
}

void order_list (struct point_list *l)
{
    uint32_t cont = 0;
    struct point_node *tmp = l->list_nodes;
    while (tmp != NULL)
    {
        tmp->id = cont;

        cont++;
        tmp = tmp->next;
    }
}

struct point_node* new_node (uint32_t id, const double pos[])
{
    struct point_node *tmp = (struct point_node*)malloc(sizeof(struct point_node));
    struct point *p = new_point(pos);
    tmp->id = id;
    tmp->value = p;
    tmp->next = NULL;
    return tmp;
}

void free_node (struct point_node *n)
{
    n->next = NULL;
    free(n);
}

struct point* new_point (const double pos[])
{
    struct point *p = (struct point*)malloc(sizeof(struct point));
    p->x = pos[0];
    p->y = pos[1];
    p->z = pos[2];
    return p;
}

void free_point (struct point *p)
{
    free(p);
    p = NULL;
}