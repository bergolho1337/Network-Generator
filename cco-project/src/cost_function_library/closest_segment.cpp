#include "closest_segment.h"

// OK
SET_COST_FUNCTION (closest_segment)
{
    struct segment_node *best = NULL;
    double minimum_length = __DBL_MAX__;

    // Get a reference to the local optimization function
    bool using_local_optimization = the_network->using_local_optimization;
    set_local_optimization_function_fn *local_optimization_fn = local_opt_config->function;

    if (using_local_optimization)
    {
        fprintf(stderr,"[closest_segment] ERROR! Local optimization is not available for the 'closest_segment'\n");
        exit(EXIT_FAILURE);
    }

    for (uint32_t i = 0; i < feasible_segments.size(); i++)
    {
        struct segment_node *iconn = feasible_segments[i];

        // First we set the bifurcation point on the middle of the segment
        struct segment_node *inew = build_segment(the_network,local_opt_config,iconn->id,new_pos);

        struct segment_node *ibiff = iconn->value->parent;

        // Calculate the length of the new segment
        double middle_pos[3];
        calc_middle_point_segment(iconn,middle_pos);

        double length = euclidean_norm(new_pos[0],new_pos[1],new_pos[2],\
                                    middle_pos[0],middle_pos[1],middle_pos[2]);

        if (length < minimum_length)
        {
            minimum_length = length;
            best = iconn;

            printf("[cost_function] Best segment = %d -- Length = %g\n",\
                                best->id,\
                                minimum_length);
        }

        restore_state_tree(the_network,iconn);
    }

    return best;
}

// OK
SET_COST_FUNCTION(closest_segment_with_length_restriction)
{

    struct segment_node *best = NULL;
    double minimum_length = __DBL_MAX__;

    // Get the length limit from the user input
    double length_limit;
    get_parameter_value_from_map(config->params,"length_limit",&length_limit);

    // Pass through the list of feasible segments
    struct segment_node *tmp;
    for (unsigned int i = 0; i < feasible_segments.size(); i++)
    {
        tmp = feasible_segments[i];
        
        double middle_pos[3];
        calc_middle_point_segment(tmp,middle_pos);

        double dist = euclidean_norm(new_pos[0],new_pos[1],new_pos[2],\
                                    middle_pos[0],middle_pos[1],middle_pos[2]);
                                    
        if (dist < minimum_length && dist >= length_limit)
        {
            minimum_length = dist;
            best = tmp;

            printf("[cost_function] Best segment = %d -- Length = %g\n",\
                                best->id,\
                                minimum_length);
        }
    }
    
    return best;
}

// OK
SET_COST_FUNCTION(closest_segment_with_angle_restriction)
{

    struct segment_node *best = NULL;
    double minimum_length = __DBL_MAX__;

    // Get parameters for the cost function
    double min_degrees_limit;
    get_parameter_value_from_map(config->params,"min_degrees_limit",&min_degrees_limit);
    double max_degrees_limit;
    get_parameter_value_from_map(config->params,"max_degrees_limit",&max_degrees_limit);

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
                                    
        if (dist < minimum_length && degrees > min_degrees_limit && degrees < max_degrees_limit)
        {
            minimum_length = dist;
            best = tmp;

            printf("[cost_function] Best segment = %d -- Length = %g -- Degrees = %g\n",\
                                best->id,\
                                minimum_length,\
                                degrees);
        }
    }
    
    return best;
}