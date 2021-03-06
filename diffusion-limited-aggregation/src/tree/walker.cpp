#include "walker.h"

struct walker* new_walker (struct user_options *the_options)
{
    struct walker *result = (struct walker*)malloc(sizeof(struct walker));

    // Get the reference to the library respawn function
    set_walker_respawn_function_fn *respawn_function_ptr = the_options->walker_config->respawn_function; 

    // Call the respawn library function
    respawn_function_ptr(the_options->walker_config,result->pos);

    result->stuck = false;
    result->stuck_counter = 0;
    
    double walker_radius = 6.0;
    get_parameter_value_from_map(the_options->walker_config->params,"walker_radius",&walker_radius);
    result->radius = walker_radius;

    return result;
}

struct walker* new_walker (struct user_options *the_options, struct walker_list *l, const uint32_t iter)
{
    //printf("Current = %u || Max = %u\n",l->num_nodes,the_options->max_num_walkers);
    struct walker *result = (struct walker*)malloc(sizeof(struct walker));

    // Get the reference to the library respawn function
    set_walker_respawn_function_fn *respawn_function_ptr = the_options->walker_config->respawn_function; 

    // Call the respawn library function
    respawn_function_ptr(the_options->walker_config,result->pos);

    // Generate a new walker position
    uint32_t counter = 0;
    bool point_is_not_ok = false;
    do
    {

        // WARNING
        if (counter > MAX_RESPAWN_TRIES)
        {
            uint32_t cur_number_walkers = l->num_nodes;
            uint32_t total_number_walkers = the_options->max_num_walkers;
            printf("[!] Warning! Reducing the number of walker to %u from %u\n",cur_number_walkers,total_number_walkers);

            the_options->max_num_walkers = cur_number_walkers;
            return NULL;
        }

        // Call the respawn library function
        respawn_function_ptr(the_options->walker_config,result->pos);

        // Check if the generated point is not part of the network
        struct walker_node *tmp = l->list_nodes;
        while (tmp != NULL)
        {
            double d = calculate_distance(result->pos,tmp->value->pos);
            if (d < EPSILON)
                point_is_not_ok = true;

            tmp = tmp->next;
        }

        counter++;
        
    }while(point_is_not_ok);

    result->stuck = false;
    result->stuck_counter = 0;
    
    double walker_radius = 6.0;
    get_parameter_value_from_map(the_options->walker_config->params,"walker_radius",&walker_radius);
    result->radius = walker_radius;
    
    return result;
}

struct walker* new_walker (const double x, const double y, const double z, const double walker_radius)
{
    struct walker *result = (struct walker*)malloc(sizeof(struct walker));

    result->pos[0] = x;
    result->pos[1] = y;
    result->pos[2] = z;

    result->stuck = false;
    result->stuck_counter = 0;
    result->radius = walker_radius;

    return result;
}

void free_walker (struct walker *the_walker)
{
    free(the_walker);
}

uint32_t is_stuck (struct walker_list *the_tree, struct walker *the_other)
{   
    struct walker_node *tmp = the_tree->list_nodes;
    double min_dist = __DBL_MAX__;
    uint32_t min_index = the_tree->num_nodes+1;
    
    // Pass through all the nodes already in the tree
    while (tmp != NULL)
    {
        struct walker *cur_node = tmp->value;
        double node_radius = cur_node->radius;
        uint32_t stuck_counter = cur_node->stuck_counter;

        // Calculate the distance between the center of the current node and the walker
        double d = calculate_distance(the_other->pos,cur_node->pos);
        // Calculate the euclidean distance between the center of the current node and the walker
        double norm = calculate_euclidean_norm(cur_node->pos[0],cur_node->pos[1],cur_node->pos[2],\
                                            the_other->pos[0],the_other->pos[1],the_other->pos[2]);

        // [RESTRICTION ZONE]
        bool is_inside = (d < node_radius * node_radius * 4.0);
        bool is_not_trifurcation = (stuck_counter < 2);

        if (is_inside && is_not_trifurcation)
        {
            // Save the node which has the minimum distance to the current walker
            if (norm < min_dist)
            {
                min_dist = norm;
                min_index = tmp->id;
            }
        }
        tmp = tmp->next;
    }

    // If the minimum distance is zero do not stuck the walker to the tree
    if (min_dist == 0.0)
        return the_tree->num_nodes+1;
    else
        return min_index;
}

void print_walker (struct walker *the_walker)
{
    printf("%g %g %g\n",the_walker->pos[0],the_walker->pos[1],the_walker->pos[2]);
}

struct walker_list* new_walker_list ()
{
    struct walker_list *s = (struct walker_list*)malloc(sizeof(struct walker_list));
    s->num_nodes = 0;
    s->list_nodes = NULL;
    return s;
}

void free_walker_list (struct walker_list *s)
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

struct walker_node* insert_walker_node (struct walker_list *l, struct walker *walker)
{
    struct walker_node *node;
    // List is empty
    if (!l->list_nodes)
    {
        node = new_walker_node(l->num_nodes,walker);
        l->list_nodes = node;
    }
    // List has some values
    else
    {
        // Pass through all the list
        struct walker_node *tmp = l->list_nodes;
        while (tmp->next != NULL)
        {
            tmp = tmp->next;
        }

        node = new_walker_node(l->num_nodes,walker);
        tmp->next = node;
    }
    l->num_nodes++;
    return node;
}

void delete_node (struct walker_list *l, const uint32_t index)
{
    struct walker_node *aux1, *aux2;
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
    free_walker(aux1->value);
    free(aux1);
    l->num_nodes--;
}

struct walker_node* search_walker_node (struct walker_list *l, const uint32_t index)
{
    struct walker_node *tmp = l->list_nodes;
    while (tmp != NULL)
    {
        if (tmp->id == index)
            return tmp;
        tmp = tmp->next;
    }
    fprintf(stderr,"[-] ERROR! walker node %u not found!\n",index);
    return NULL;
}

void order_list (struct walker_list *l)
{
    uint32_t cont = 0;
    struct walker_node *tmp = l->list_nodes;
    while (tmp != NULL)
    {
        tmp->id = cont;

        cont++;
        tmp = tmp->next;
    }
}

bool is_empty (struct walker_list *l)
{
    return (l->list_nodes == NULL) ? true : false;
}

bool is_duplicate (const double pos[], struct walker_list *l)
{
    struct walker_node *tmp = l->list_nodes;
    while (tmp != NULL)
    {
        double d = calculate_distance(tmp->value->pos,pos);
        if (d < EPSILON)
            return true;

        tmp = tmp->next;
    }
    return false;
}

struct walker_node* new_walker_node (uint32_t id, struct walker *s)
{
    struct walker_node *tmp = (struct walker_node*)malloc(sizeof(struct walker_node));
    tmp->id = id;
    tmp->value = s;
    tmp->next = NULL;
    return tmp;
}

void print_list (struct walker_list *l)
{
    if (l->num_nodes == 0)
    {
        printf("[!] The list is empty!\n");
        return;
    }

    struct walker_node *tmp = l->list_nodes;
    printf("Number of walker_node in list = %u\n",l->num_nodes);
    while (tmp != NULL)
    {
        printf("Walker %u - (%g,%g,%g) || Links = %u\n",tmp->id,tmp->value->pos[0],tmp->value->pos[1],tmp->value->pos[2],tmp->value->stuck_counter);
        tmp = tmp->next;
    }
}

void write_list (struct walker_list *l, const uint32_t iter)
{
    uint32_t num_points = l->num_nodes;

    char filename[50];
    sprintf(filename,"output/walker/walker-%u.vtk",iter);
    FILE *file = fopen(filename,"w+");

    fprintf(file,"# vtk DataFile Version 3.0\n");
    fprintf(file,"Walker\n");
    fprintf(file,"ASCII\n");
    fprintf(file,"DATASET POLYDATA\n");
    fprintf(file,"POINTS %u float\n",num_points);
    
    // Write the points
    struct walker_node *tmp = l->list_nodes;
    while (tmp != NULL)
    {
        fprintf(file,"%g %g %g\n",tmp->value->pos[0],\
                                tmp->value->pos[1],\
                                tmp->value->pos[2]);
        tmp = tmp->next;
    }
    fprintf(file,"VERTICES %lu %lu\n",num_points,num_points*2);    
    for (uint32_t i = 0; i < num_points; i++)
        fprintf(file,"1 %u\n",i);
    
    fclose(file);
}