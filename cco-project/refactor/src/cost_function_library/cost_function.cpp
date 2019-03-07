#include "cost_function.h"

SET_COST_FUNCTION(funcao1)
{
    printf("Funcao 1: %s\n",config->name);

    double value;
    get_parameter_value_from_map(config->params,"param1",&value);
    printf("Value = %g\n",value);
    
    return NULL;
}

SET_COST_FUNCTION(funcao2)
{
    printf("Funcao 2: %s\n",config->name);

    return NULL;
}

SET_COST_FUNCTION(funcao3)
{
    printf("Funcao 3: %s\n",config->name);

    return NULL;
}

SET_COST_FUNCTION(funcao4)
{
    printf("Funcao 4: %s\n",config->name);

    return NULL;
}

SET_COST_FUNCTION(closest_segment)
{
    
    struct segment_node *closest = NULL;
    double closest_dist = __DBL_MAX__;

    // Pass through the list of feasible segments
    struct segment_node *tmp;
    for (unsigned int i = 0; i < feasible_segments.size(); i++)
    {
        tmp = feasible_segments[i];
        
        double middle_pos[3];
        calc_middle_point_segment(tmp,middle_pos);

        double dist = euclidean_norm(new_pos[0],new_pos[1],new_pos[2],\
                                    middle_pos[0],middle_pos[1],middle_pos[2]);
        if (dist < closest_dist)
        {
            closest_dist = dist;
            closest = tmp;
        }
    }
    
    return closest;
}

SET_COST_FUNCTION(closest_segment_with_limit_size)
{

    struct segment_node *closest = NULL;
    double closest_dist = __DBL_MAX__;

    // Get parameters for the cost function
    double size_limit;
    get_parameter_value_from_map(config->params,"size_limit",&size_limit);

    // Pass through the list of feasible segments
    struct segment_node *tmp;
    for (unsigned int i = 0; i < feasible_segments.size(); i++)
    {
        tmp = feasible_segments[i];
        
        double middle_pos[3];
        calc_middle_point_segment(tmp,middle_pos);

        double dist = euclidean_norm(new_pos[0],new_pos[1],new_pos[2],\
                                    middle_pos[0],middle_pos[1],middle_pos[2]);
                                    
        if (dist < closest_dist && dist >= size_limit)
        {
            closest_dist = dist;
            closest = tmp;
        }
    }
    
    return closest;
}

SET_COST_FUNCTION(closest_segment_with_angle_restriction)
{

    struct segment_node *closest = NULL;
    double closest_dist = __DBL_MAX__;

    // Get parameters for the cost function
    double degrees_limit;
    get_parameter_value_from_map(config->params,"degrees",&degrees_limit);

    // Pass through the list of feasible segments
    struct segment_node *tmp;
    for (unsigned int i = 0; i < feasible_segments.size(); i++)
    {
        tmp = feasible_segments[i];
        
        struct point *src = tmp->value->src->value;
        struct point *dest = tmp->value->dest->value;

        double middle_pos[3];
        calc_middle_point_segment(tmp,middle_pos);

        double u[3], v[3];
        build_unitary_vector(u,middle_pos[0],middle_pos[1],middle_pos[2],\
                              dest->x,dest->y,dest->z);
        build_unitary_vector(v,middle_pos[0],middle_pos[1],middle_pos[2],\
                              new_pos[0],new_pos[1],new_pos[2]);

        double dist = euclidean_norm(new_pos[0],new_pos[1],new_pos[2],\
                                    middle_pos[0],middle_pos[1],middle_pos[2]);

        double degrees = calc_angle_between_vectors(u,v);
                                    
        if (dist < closest_dist && degrees >= degrees_limit)
        {
            closest_dist = dist;
            closest = tmp;
        }
    }
    
    return closest;
}

SET_COST_FUNCTION (minimize_tree_volume)
{
    struct segment_node *best = NULL;
    double minimum_volume = __DBL_MAX__;

    for (uint32_t i = 0; i < feasible_segments.size(); i++)
    {
        struct segment_node *iconn = feasible_segments[i];

        //printf("\nBefore\n");
        //print_list(the_network->segment_list);
        //printf("\n");

        build_segment(the_network,iconn->id,new_pos);

        double volume = calc_tree_volume(the_network);
        if (volume < minimum_volume)
        {
            minimum_volume = volume;
            best = iconn;

            printf("[cost_function] Best segment = %d -- Volume = %g\n",best->id,minimum_volume);
        }

        //printf("\nProcessing\n");
        //print_list(the_network->segment_list);
        //printf("\n");

        restore_state_tree(the_network,iconn);

        //printf("\nAfter\n");
        //print_list(the_network->segment_list);
        //printf("\n");
    }

    return best;
}