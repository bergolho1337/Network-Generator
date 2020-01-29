#include "segment.h"

struct segment* new_segment (const uint32_t src, const uint32_t dest)
{
    struct segment *result = (struct segment*)malloc(sizeof(struct segment));

    result->src = src;
    result->dest = dest;

    return result;
}

void free_segment (struct segment *the_segment)
{
    free(the_segment);
}

struct segment_list* new_segment_list ()
{
    struct segment_list *s = (struct segment_list*)malloc(sizeof(struct segment_list));
    s->num_nodes = 0;
    s->list_nodes = NULL;
    return s;
}

void free_segment_list (struct segment_list *s)
{
    uint32_t cont = 0;
    while (!is_empty(s))
    {
        delete_node(s,cont);
        cont++;
    }
    //print_list(l);
    
    free(s);
    s = NULL;
}

struct segment_node* insert_segment_node (struct segment_list *l, struct segment *segment)
{
    struct segment_node *node;
    // List is empty
    if (!l->list_nodes)
    {
        node = new_segment_node(l->num_nodes,segment);
        l->list_nodes = node;
    }
    // List has some values
    else
    {
        // Pass through all the list
        struct segment_node *tmp = l->list_nodes;
        while (tmp->next != NULL)
        {
            tmp = tmp->next;
        }

        node = new_segment_node(l->num_nodes,segment);
        tmp->next = node;
    }
    l->num_nodes++;
    return node;
}

void delete_node (struct segment_list *l, const uint32_t index)
{
    struct segment_node *aux1, *aux2;
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
    free_segment(aux1->value);
    free(aux1);
    l->num_nodes--;
}

struct segment_node* search_segment_node (struct segment_list *l, const uint32_t index)
{
    struct segment_node *tmp = l->list_nodes;
    while (tmp != NULL)
    {
        if (tmp->id == index)
            return tmp;
        tmp = tmp->next;
    }
    fprintf(stderr,"[-] ERROR! Segment node %u not found!\n",index);
    return NULL;
}

bool is_empty (struct segment_list *l)
{
    return (l->list_nodes == NULL) ? true : false;
}

struct segment_node* new_segment_node (uint32_t id, struct segment *s)
{
    struct segment_node *tmp = (struct segment_node*)malloc(sizeof(struct segment_node));
    tmp->id = id;
    tmp->value = s;
    tmp->next = NULL;
    return tmp;
}

void print_list (struct segment_list *l)
{
    if (l->num_nodes == 0)
    {
        printf("[!] The list is empty!\n");
        return;
    }

    struct segment_node *tmp = l->list_nodes;
    printf("Number of segment_node in list = %u\n",l->num_nodes);
    while (tmp != NULL)
    {
        printf("Segment %d (%d,%d)\n",tmp->id,tmp->value->src,tmp->value->dest);
        tmp = tmp->next;
    }
}