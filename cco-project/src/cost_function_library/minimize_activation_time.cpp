#include "minimize_activation_time.h"

SET_COST_FUNCTION (minimize_tree_activation_time)
{
    struct segment_node *best = NULL;
    double minimum_at = __DBL_MAX__;

    double Q_perf = the_network->Q_perf;
    double p_perf = the_network->p_perf;
    double p_term = the_network->p_term;
    double delta_p = p_perf - p_term;

    // Get a reference to the local optimization function
    bool using_local_optimization = the_network->using_local_optimization;
    set_local_optimization_function_fn *local_optimization_fn = local_opt_config->function;

    // Get cost function parameters
    bool ret;
    double c;
    ret = get_parameter_value_from_map(config->params,"c",&c);
    double cm;
    ret = get_parameter_value_from_map(config->params,"cm",&cm);
    double rc;
    ret = get_parameter_value_from_map(config->params,"rc",&rc);
    double rm;
    ret = get_parameter_value_from_map(config->params,"rm",&rm);
    double deviation_limit = __DBL_MAX__;
    ret = get_parameter_value_from_map(config->params,"deviation_limit",&deviation_limit);
    double min_angle_limit = 0.0;
    ret = get_parameter_value_from_map(config->params,"min_angle_limit",&min_angle_limit);
    double max_angle_limit = 360.0;
    ret = get_parameter_value_from_map(config->params,"max_angle_limit",&max_angle_limit);

    //printf("Feasible segment %u\n",feasible_segments.size());

    for (uint32_t i = 0; i < feasible_segments.size(); i++)
    {

        struct segment_node *iconn = feasible_segments[i];

        struct segment_node *inew = build_segment(the_network,local_opt_config,iconn->id,new_pos);

        struct segment_node *ibiff = iconn->value->parent;

        //  WITH LOCAL OPTIMIZATION
        if (using_local_optimization)
        {

            // Save the original position of the bifurcation
            double ori_pos[3];
            save_original_bifurcation_position(iconn,ori_pos);

            // Initialize the best position as middle point of the segment
            double best_pos[3];
            initialize_best_position_as_middle_point(best_pos,ori_pos);

            //printf("[cost_function] Original bifurcation position = (%g,%g,%g)\n",\
                                                        best_pos[0],\
                                                        best_pos[1],\
                                                        best_pos[2]);

            // 1) Check the cost function of the first configuration
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

            double angle = calc_angle_between_vectors(u,v);

            if (at < minimum_at && \
            !has_deviation(the_network->segment_list,inew,at,deviation_limit,c,cm,rc,rm) && \
            angle >= min_angle_limit && angle <= max_angle_limit)
            {
                minimum_at = at;
                best = iconn;

                // The best position of the best segment will be stored inside the
                // 'local_optimization' structure
                local_opt_config->best_pos[0] = best_pos[0];
                local_opt_config->best_pos[1] = best_pos[1];
                local_opt_config->best_pos[2] = best_pos[2];

                printf("[cost_function] Best segment = %d -- Activation time = %g ms -- Best position = (%g,%g,%g)\n",\
                                best->id,\
                                minimum_at,\
                                local_opt_config->best_pos[0],\
                                local_opt_config->best_pos[1],\
                                local_opt_config->best_pos[2]);

            }

            // 2) Now, call the local optimization function and fill the 'test_positions' array
            std::vector<struct point> test_positions;
            local_optimization_fn(iconn,ibiff,inew,test_positions);

            for (uint32_t j = 0; j < test_positions.size(); j++)
            {
                //printf("Test position = (%g,%g,%g)\n",test_positions[j].x,test_positions[j].y,test_positions[j].z);

                // Change the position of the bifurcation point
                double new_biff_pos[3];
                new_biff_pos[0] = test_positions[j].x;
                new_biff_pos[1] = test_positions[j].y;
                new_biff_pos[2] = test_positions[j].z;
                move_bifurcation_location(iconn,ibiff,inew,new_biff_pos);
                rescale_tree(ibiff,iconn,inew,Q_perf,delta_p,the_network->num_terminals);
                recalculate_radius(the_network);

                // Evaluate cost function for the current configuration
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

                double angle = calc_angle_between_vectors(u,v);

                if (at < minimum_at && \
            !has_deviation(the_network->segment_list,inew,at,deviation_limit,c,cm,rc,rm) && \
            angle >= min_angle_limit && angle <= max_angle_limit)
                {
                    minimum_at = at;
                    best = iconn;

                    // The best position of the best segment will be stored inside the 'local_optimization' structure
                    local_opt_config->best_pos[0] = new_biff_pos[0];
                    local_opt_config->best_pos[1] = new_biff_pos[1];
                    local_opt_config->best_pos[2] = new_biff_pos[2];

                    printf("[cost_function] Best segment = %d -- Activation time = %g ms -- Best position = (%g,%g,%g)\n",\
                                    best->id,\
                                    minimum_at,\
                                    local_opt_config->best_pos[0],\
                                    local_opt_config->best_pos[1],\
                                    local_opt_config->best_pos[2]);
                }

            }

            // Move the bifurcation to the original position
            move_bifurcation_location(iconn,ibiff,inew,ori_pos);

            restore_state_tree(the_network,iconn);

            // Clear the array for the next segment
            test_positions.clear();
 
        }
        // NO LOCAL OPTIMIZATION
        else
        {
            // Evaluate the cost function
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

            double angle = calc_angle_between_vectors(u,v);

            if (at < minimum_at && \
            !has_deviation(the_network->segment_list,inew,at,deviation_limit,c,cm,rc,rm) && \
            angle >= min_angle_limit && angle <= max_angle_limit)
            {
                minimum_at = at;
                best = iconn;

                printf("[cost_function] Best segment = %d -- Activation time = %g ms\n",best->id,minimum_at);
            }

            restore_state_tree(the_network,iconn);
        }
    }

    if (using_local_optimization)
    {
        // Set off the 'first_call' flag
        local_opt_config->first_call = false;
    }

    // DEBUG
    //printf("[cost_function] Best AT = %g ms\n",minimum_at);
    //print_terminal_activation_time(the_network,c,cm,rc,rm);

    return best;
}

// TODO: Need to implement ...
/*
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
    printf("[cost_function] Best AT = %g ms\n",minimum_at);
    //print_terminal_activation_time(the_network,c,cm,rc,rm);

    return best;
}

SET_COST_FUNCTION (minimize_tree_activation_time_with_angle_restriction_and_level_restriction)
{
    struct segment_node *best = NULL;
    double minimum_eval = __DBL_MAX__;

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
    double degrees_min;
    get_parameter_value_from_map(config->params,"degrees_min",&degrees_min);
    double degrees_max;
    get_parameter_value_from_map(config->params,"degrees_max",&degrees_max);

    for (uint32_t i = 0; i < feasible_segments.size(); i++)
    {
        struct segment_node *iconn = feasible_segments[i];

        struct segment_node *inew = build_segment(the_network,iconn->id,new_pos);

        // Calculate activation time
        double at = calc_terminal_activation_time(inew,c,cm,rc,rm);

        // Convert from miliseconds to microseconds
        double new_at = at * MS_TO_US;

        // Call the segment level function
        double eval = calc_segment_activation_time_using_level(new_at,iconn);

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

        if (eval < minimum_eval && \
            !has_deviation(the_network->segment_list,inew,at,deviation_limit,c,cm,rc,rm) && \
            degrees >= degrees_min && degrees <= degrees_max)
        {
            best = iconn;
            minimum_eval = eval;
        }

        restore_state_tree(the_network,iconn);
    }

    // DEBUG
    printf("[cost_function] Best eval = %g us\n",minimum_eval);
    //print_terminal_activation_time(the_network,c,cm,rc,rm);

    return best;
}
*/
