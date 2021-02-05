#include "linked-list.h"

struct linked_list* new_linked_list ()
{
    struct linked_list *l = (struct linked_list*)malloc(sizeof(struct linked_list));
    l->num_nodes = 0;
    l->list_nodes = NULL;
    return l;
}

void free_linked_list (struct linked_list *l)
{
    while (!is_empty(l))
    {
        NODE_DATATYPE key = l->list_nodes->value;
        delete_node(l,key);
    }
    //print_list(l);
    
    free(l);
    l = NULL;
}

void insert_node (struct linked_list *l, const NODE_DATATYPE key)
{
    // List is empty
    if (!l->list_nodes)
    {
        struct node *node = new_node(l->num_nodes,key);
        l->list_nodes = node;
    }
    // List has some values
    else
    {
        // Pass through all the list
        struct node *tmp = l->list_nodes;
        while (tmp->next != NULL)
        {
            tmp = tmp->next;
        }

        struct node *node = new_node(l->num_nodes,key);
        tmp->next = node;
    }
    l->num_nodes++;
}

void delete_node (struct linked_list *l, const NODE_DATATYPE key)
{
    struct node *aux1, *aux2;
    aux1 = l->list_nodes;
    aux2 = NULL;

    if (l->num_nodes == 0)
    {
        fprintf(stderr,"[-] ERROR! The list is empty!\n");
        exit(EXIT_FAILURE);
    }

    while (aux1 != NULL && aux1->value != key)
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
    free(aux1);
    l->num_nodes--;
}

struct node* search_node (struct linked_list *l, const NODE_DATATYPE key)
{
    struct node *tmp = l->list_nodes;
    while (tmp != NULL)
    {
        if (tmp->value == key)
        {
            printf("[+] Found the node || %u %lf ||\n",tmp->id,tmp->value);
            return tmp;
        }
            
        tmp = tmp->next;
    }
    fprintf(stderr,"[-] Key %lf not found!\n",key);
    return NULL;
}

bool is_empty (struct linked_list *l)
{
    return (l->list_nodes == NULL) ? true : false;
}

void print_list (struct linked_list *l)
{
    if (l->num_nodes == 0)
    {
        printf("[!] The list is empty!\n");
        return;
    }

    struct node *tmp = l->list_nodes;
    printf("Number of node in list = %u\n",l->num_nodes);
    while (tmp != NULL)
    {
        printf("Node %d -- Value = %g\n",tmp->id,tmp->value);
        tmp = tmp->next;
    }
}

void order_list (struct linked_list *l)
{
    uint32_t cont = 0;
    struct node *tmp = l->list_nodes;
    while (tmp != NULL)
    {
        tmp->id = cont;

        cont++;
        tmp = tmp->next;
    }
}

struct node* new_node (uint32_t id, const NODE_DATATYPE key)
{
    struct node *tmp = (struct node*)malloc(sizeof(struct node));
    tmp->id = id;
    tmp->value = key;
    tmp->next = NULL;
    return tmp;
}

void free_node (struct node *n)
{
    n->next = NULL;
    free(n);
}