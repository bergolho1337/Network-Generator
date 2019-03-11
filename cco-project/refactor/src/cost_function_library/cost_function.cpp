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

    double epsilon_lim;
    get_parameter_value_from_map(config->params,"epsilon_lim",&epsilon_lim);

    for (uint32_t i = 0; i < feasible_segments.size(); i++)
    {
        struct segment_node *iconn = feasible_segments[i];

        struct segment_node *inew = build_segment(the_network,iconn->id,new_pos);

        double epsilon_rad = calc_assymetric_ratio(iconn,inew);

        double volume = calc_tree_volume(the_network);
        if (volume < minimum_volume && epsilon_rad > epsilon_lim)
        {
            minimum_volume = volume;
            best = iconn;

            printf("[cost_function] Best segment = %d -- Volume = %g\n",best->id,minimum_volume);
        }

        restore_state_tree(the_network,iconn);

    }

    return best;
}

SET_COST_FUNCTION (maximize_tree_volume)
{
    struct segment_node *best = NULL;
    double maximum_volume = __DBL_MIN__;

    for (uint32_t i = 0; i < feasible_segments.size(); i++)
    {
        struct segment_node *iconn = feasible_segments[i];

        build_segment(the_network,iconn->id,new_pos);

        double volume = calc_tree_volume(the_network);
        if (volume > maximum_volume)
        {
            maximum_volume = volume;
            best = iconn;

            printf("[cost_function] Best segment = %d -- Volume = %g\n",best->id,maximum_volume);
        }

        restore_state_tree(the_network,iconn);

    }

    return best;
}

SET_COST_FUNCTION (minimize_tree_activation_time)
{
    struct segment_node *best = NULL;
    double minimum_at = __DBL_MAX__;

    // Get cost function parameters
    double c;
    get_parameter_value_from_map(config->params,"c",&c);
    double cm;
    get_parameter_value_from_map(config->params,"cm",&cm);
    double rc;
    get_parameter_value_from_map(config->params,"rc",&rc);
    double rm;
    get_parameter_value_from_map(config->params,"rm",&rm);
    double deviation_limit;
    get_parameter_value_from_map(config->params,"deviation_limit",&deviation_limit);

    for (uint32_t i = 0; i < feasible_segments.size(); i++)
    {
        struct segment_node *iconn = feasible_segments[i];

        struct segment_node *inew = build_segment(the_network,iconn->id,new_pos);

        double at = calc_terminal_activation_time(inew,c,cm,rc,rm);
        
        if (at < minimum_at && \
            !has_deviation(the_network->segment_list,inew,at,deviation_limit,c,cm,rc,rm))
        {
            best = iconn;
            minimum_at = at;
        }

        restore_state_tree(the_network,iconn);
    }

    // DEBUG
    printf("[cost_function] Best AT = %g\n",minimum_at);
    print_terminal_activation_time(the_network,c,cm,rc,rm);

    return best;
}

SET_COST_FUNCTION (minimize_tree_activation_time_with_angle_restriction)
{
    struct segment_node *best = NULL;
    double minimum_at = __DBL_MAX__;

    // Get cost function parameters
    double c;
    get_parameter_value_from_map(config->params,"c",&c);
    double cm;
    get_parameter_value_from_map(config->params,"cm",&cm);
    double rc;
    get_parameter_value_from_map(config->params,"rc",&rc);
    double rm;
    get_parameter_value_from_map(config->params,"rm",&rm);
    double deviation_limit;
    get_parameter_value_from_map(config->params,"deviation_limit",&deviation_limit);
    double degrees_limit;
    get_parameter_value_from_map(config->params,"degrees",&degrees_limit);

    for (uint32_t i = 0; i < feasible_segments.size(); i++)
    {
        struct segment_node *iconn = feasible_segments[i];

        struct segment_node *inew = build_segment(the_network,iconn->id,new_pos);

        // Calculate activation time
        double at = calc_terminal_activation_time(inew,c,cm,rc,rm);

        // Calculate bifurcation angle
        struct point *src = iconn->value->src->value;
        struct point *dest = iconn->value->dest->value;

        double middle_pos[3];
        calc_middle_point_segment(iconn,middle_pos);

        double u[3], v[3];
        build_unitary_vector(u,middle_pos[0],middle_pos[1],middle_pos[2],\
                              dest->x,dest->y,dest->z);
        build_unitary_vector(v,middle_pos[0],middle_pos[1],middle_pos[2],\
                              new_pos[0],new_pos[1],new_pos[2]);

        double degrees = calc_angle_between_vectors(u,v);

        if (at < minimum_at && \
            !has_deviation(the_network->segment_list,inew,at,deviation_limit,c,cm,rc,rm) && \
            degrees >= degrees_limit)
        {
            best = iconn;
            minimum_at = at;
        }

        restore_state_tree(the_network,iconn);
    }

    // DEBUG
    printf("[cost_function] Best AT = %g\n",minimum_at);
    print_terminal_activation_time(the_network,c,cm,rc,rm);

    return best;
}