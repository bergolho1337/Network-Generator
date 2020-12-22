#include "face-list.h"

struct face_list* new_face_list ()
{
    struct face_list *l = (struct face_list*)malloc(sizeof(struct face_list));
    l->num_nodes = 0;
    l->list_nodes = NULL;
    return l;
}

void free_face_list (struct face_list *l)
{
    uint32_t cont = 0;
    while (!is_empty(l))
    {
        delete_node(l,0);
        cont++;
    }
    //print_list(l);
    
    free(l);
    l = NULL;
}

struct face_node* insert_point (struct face_list *l,\
                            const double pos1[], const double pos2[], const double pos3[], const double n[])
{
    struct face_node *node;
    // List is empty
    if (!l->list_nodes)
    {
        node = new_node(l->num_nodes,pos1,pos2,pos3,n);
        l->list_nodes = node;
    }
    // List has some values
    else
    {
        // Pass through all the list
        struct face_node *tmp = l->list_nodes;
        while (tmp->next != NULL)
        {
            tmp = tmp->next;
        }

        node = new_node(l->num_nodes,pos1,pos2,pos3,n);
        tmp->next = node;
    }
    l->num_nodes++;
    return node;
}

void delete_node (struct face_list *l, const uint32_t index)
{
    struct face_node *aux1, *aux2;
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
    free_face(aux1->value);
    free(aux1);
    l->num_nodes--;

    if (l)
        order_list(l);
}

struct face_node* search_node (struct face_list *l, const uint32_t index)
{
    struct face_node *tmp = l->list_nodes;
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

bool is_empty (struct face_list *l)
{
    return (l->list_nodes == NULL) ? true : false;
}

void print_list (struct face_list *l)
{
    if (l->num_nodes == 0)
    {
        printf("[!] The list is empty!\n");
        return;
    }

    struct face_node *tmp = l->list_nodes;
    printf("Number of node in list = %u\n",l->num_nodes);
    while (tmp != NULL)
    {
        printf("Node %d -- Vertex_1(%g,%g,%g) || Vertex_2(%g,%g,%g) || Vertex_3(%g,%g,%g) || Normal(%g,%g,%g)\n",tmp->id,\
                                        tmp->value->x1,tmp->value->y1,tmp->value->z1,\
                                        tmp->value->x2,tmp->value->y2,tmp->value->z2,\
                                        tmp->value->x3,tmp->value->y3,tmp->value->z3,\
                                        tmp->value->nx,tmp->value->ny,tmp->value->nz);
        tmp = tmp->next;
    }
}

void write_list (struct face_list *l, FILE *log_file)
{
    if (l->num_nodes == 0)
    {
        fprintf(log_file,"[!] The list is empty!\n");
        return;
    }

    struct face_node *tmp = l->list_nodes;
    fprintf(log_file,"Number of node in list = %u\n",l->num_nodes);
    while (tmp != NULL)
    {
        fprintf(log_file,"Node %d -- Vertex_1(%g,%g,%g) || Vertex_2(%g,%g,%g) || Vertex_3(%g,%g,%g) || Normal(%g,%g,%g)\n",tmp->id,\
                                        tmp->value->x1,tmp->value->y1,tmp->value->z1,\
                                        tmp->value->x2,tmp->value->y2,tmp->value->z2,\
                                        tmp->value->x3,tmp->value->y3,tmp->value->z3,\
                                        tmp->value->nx,tmp->value->ny,tmp->value->nz);
        tmp = tmp->next;
    }
}

void order_list (struct face_list *l)
{
    uint32_t cont = 0;
    struct face_node *tmp = l->list_nodes;
    while (tmp != NULL)
    {
        tmp->id = cont;

        cont++;
        tmp = tmp->next;
    }
}

struct face_node* new_node (uint32_t id, const double pos1[], const double pos2[], const double pos3[], const double n[])
{
    struct face_node *tmp = (struct face_node*)malloc(sizeof(struct face_node));
    struct face *f = new_face(pos1,pos2,pos3,n);
    tmp->id = id;
    tmp->value = f;
    tmp->next = NULL;
    return tmp;
}

void free_node (struct face_node *n)
{
    n->next = NULL;
    free(n);
}

struct face* new_face (const double pos1[], const double pos2[], const double pos3[], const double n[])
{
    struct face *f = (struct face*)malloc(sizeof(struct face));
    f->x1 = pos1[0]; f->y1 = pos1[1]; f->z1 = pos1[2];
    f->x2 = pos2[0]; f->y2 = pos2[1]; f->z2 = pos2[2];
    f->x3 = pos3[0]; f->y3 = pos3[1]; f->z3 = pos3[2];
    f->nx = n[0]; f->ny = n[1]; f->nz = n[2];
    return f;
}

void free_point (struct point *p)
{
    free(p);
    p = NULL;
}