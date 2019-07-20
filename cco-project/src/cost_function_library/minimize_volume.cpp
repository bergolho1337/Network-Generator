#include "minimize_volume.h"

// TODO: Move all minimize_tree_volume functions to a separate file ...
SET_COST_FUNCTION (minimize_tree_volume_default)
{
    struct segment_node *best = NULL;
    double minimum_volume = __DBL_MAX__;

    double Q_perf = the_network->Q_perf;
    double p_perf = the_network->p_perf;
    double p_term = the_network->p_term;
    double delta_p = p_perf - p_term;

    // Get a reference to the local optimization function
    bool using_local_optimization = the_network->using_local_optimization;
    set_local_optimization_function_fn *local_optimization_fn = local_opt_config->function;

    for (uint32_t i = 0; i < feasible_segments.size(); i++)
    {
        struct segment_node *iconn = feasible_segments[i];

        // First we set the bifurcation point on the middle of the segment
        struct segment_node *inew = build_segment(the_network,local_opt_config,iconn->id,new_pos);

        struct segment_node *ibiff = iconn->value->parent;

        //  WITH LOCAL OPTIMIZATION
        if (using_local_optimization)
        {
            // Save the original position of the bifurcation
            double ori_pos[3];
            save_original_bifurcation_position(ibiff,ori_pos);

            // Initialize the best position as middle point of the segment
            double best_pos[3];
            initialize_best_position_as_middle_point(best_pos,ori_pos);

            //printf("[cost_function] Original bifurcation position = (%g,%g,%g)\n",\
                                                        best_pos[0],\
                                                        best_pos[1],\
                                                        best_pos[2]);

            // 1) Check the cost function of the first configuration
            double volume = calc_tree_volume(the_network);
            if (volume < minimum_volume)
            {
                minimum_volume = volume;
                best = iconn;

                // The best position of the best segment will be stored inside the 
                // 'local_optimization' structure
                local_opt_config->best_pos[0] = best_pos[0];
                local_opt_config->best_pos[1] = best_pos[1];
                local_opt_config->best_pos[2] = best_pos[2];

                printf("[cost_function] Best segment = %d -- Volume = %g -- Best position = (%g,%g,%g)\n",\
                                best->id,\
                                minimum_volume,\
                                local_opt_config->best_pos[0],\
                                local_opt_config->best_pos[1],\
                                local_opt_config->best_pos[2]);

            }

            // 2) Now, call the local optimization function and fill the 'test_positions' array
            std::vector<struct point> test_positions;
            local_optimization_fn(iconn,ibiff,inew,test_positions);

            for (uint32_t j = 0; j < test_positions.size(); j++)
            {
                // Change the position of the bifurcation point 
                double new_pos[3];
                new_pos[0] = test_positions[j].x;
                new_pos[1] = test_positions[j].y;
                new_pos[2] = test_positions[j].z;
                move_bifurcation_location(iconn,ibiff,inew,new_pos);

                rescale_tree(ibiff,iconn,inew,Q_perf,delta_p,the_network->num_terminals);

                recalculate_radius(the_network);

                double volume = calc_tree_volume(the_network);
                if (volume < minimum_volume)
                {
                    minimum_volume = volume;
                    best = iconn;

                    // The best position of the best segment will be stored inside the 'local_optimization' structure
                    local_opt_config->best_pos[0] = new_pos[0];
                    local_opt_config->best_pos[1] = new_pos[1];
                    local_opt_config->best_pos[2] = new_pos[2];

                    printf("[cost_function] Best segment = %d -- Volume = %g -- Best position = (%g,%g,%g)\n",\
                                    best->id,\
                                    minimum_volume,\
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
            double volume = calc_tree_volume(the_network);
            if (volume < minimum_volume)
            {
                minimum_volume = volume;
                best = iconn;

                printf("[cost_function] Best segment = %d -- Volume = %g\n",best->id,minimum_volume);
            }

            restore_state_tree(the_network,iconn);
        }
    }

    return best;
}

SET_COST_FUNCTION (minimize_tree_volume_with_angle_restriction)
{
    struct segment_node *best = NULL;
    double minimum_volume = __DBL_MAX__;

    double Q_perf = the_network->Q_perf;
    double p_perf = the_network->p_perf;
    double p_term = the_network->p_term;
    double delta_p = p_perf - p_term;

    // Get a reference to the local optimization function
    bool using_local_optimization = the_network->using_local_optimization;
    set_local_optimization_function_fn *local_optimization_fn = local_opt_config->function;

    // Get the degree limit from the user input
    double degrees_limit;
    get_parameter_value_from_map(config->params,"degrees",&degrees_limit);

    for (uint32_t i = 0; i < feasible_segments.size(); i++)
    {
        struct segment_node *iconn = feasible_segments[i];

        // First we set the bifurcation point on the middle of the segment
        struct segment_node *inew = build_segment(the_network,local_opt_config,iconn->id,new_pos);

        struct segment_node *ibiff = iconn->value->parent;

        //  WITH LOCAL OPTIMIZATION
        if (using_local_optimization)
        {
            // Save the original position of the bifurcation
            double ori_pos[3];
            save_original_bifurcation_position(ibiff,ori_pos);

            // Initialize the best position as middle point of the segment
            double best_pos[3];
            initialize_best_position_as_middle_point(best_pos,ori_pos);

            //printf("[cost_function] Original bifurcation position = (%g,%g,%g)\n",\
                                                        best_pos[0],\
                                                        best_pos[1],\
                                                        best_pos[2]);

            // 1) Check the cost function of the first configuration
            double volume = calc_tree_volume(the_network);

            // 2) Calculate bifurcation angle
            struct point *src = iconn->value->src->value;
            struct point *dest = iconn->value->dest->value;

            double middle_pos[3];
            calc_middle_point_segment(iconn,middle_pos);

            double u[3], v[3];
            build_unitary_vector(u,middle_pos[0],middle_pos[1],middle_pos[2],\
                                dest->x,dest->y,dest->z);
            build_unitary_vector(v,middle_pos[0],middle_pos[1],middle_pos[2],\
                                new_pos[0],new_pos[1],new_pos[2]);

            double dist = euclidean_norm(new_pos[0],new_pos[1],new_pos[2],\
                                        middle_pos[0],middle_pos[1],middle_pos[2]);

            double degrees = calc_angle_between_vectors(u,v);

            if (volume < minimum_volume && degrees > degrees_limit)
            {
                minimum_volume = volume;
                best = iconn;

                // The best position of the best segment will be stored inside the 
                // 'local_optimization' structure
                local_opt_config->best_pos[0] = best_pos[0];
                local_opt_config->best_pos[1] = best_pos[1];
                local_opt_config->best_pos[2] = best_pos[2];

                printf("[cost_function] Best segment = %d -- Volume = %g -- Best position = (%g,%g,%g)\n",\
                                best->id,\
                                minimum_volume,\
                                local_opt_config->best_pos[0],\
                                local_opt_config->best_pos[1],\
                                local_opt_config->best_pos[2]);

            }

            // 2) Now, call the local optimization function and fill the 'test_positions' array
            std::vector<struct point> test_positions;
            local_optimization_fn(iconn,ibiff,inew,test_positions);

            for (uint32_t j = 0; j < test_positions.size(); j++)
            {
                // Change the position of the bifurcation point 
                double new_pos[3];
                new_pos[0] = test_positions[j].x;
                new_pos[1] = test_positions[j].y;
                new_pos[2] = test_positions[j].z;
                move_bifurcation_location(iconn,ibiff,inew,new_pos);

                rescale_tree(ibiff,iconn,inew,Q_perf,delta_p,the_network->num_terminals);

                recalculate_radius(the_network);

                // Volume
                double volume = calc_tree_volume(the_network);

                // Angle
                struct point *src = iconn->value->src->value;
                struct point *dest = iconn->value->dest->value;

                double middle_pos[3];
                calc_middle_point_segment(iconn,middle_pos);

                double u[3], v[3];
                build_unitary_vector(u,middle_pos[0],middle_pos[1],middle_pos[2],\
                                    dest->x,dest->y,dest->z);
                build_unitary_vector(v,middle_pos[0],middle_pos[1],middle_pos[2],\
                                    new_pos[0],new_pos[1],new_pos[2]);

                double dist = euclidean_norm(new_pos[0],new_pos[1],new_pos[2],\
                                            middle_pos[0],middle_pos[1],middle_pos[2]);

                double degrees = calc_angle_between_vectors(u,v);

                if (volume < minimum_volume && degrees > degrees_limit)
                {
                    minimum_volume = volume;
                    best = iconn;

                    // The best position of the best segment will be stored inside the 'local_optimization' structure
                    local_opt_config->best_pos[0] = new_pos[0];
                    local_opt_config->best_pos[1] = new_pos[1];
                    local_opt_config->best_pos[2] = new_pos[2];

                    printf("[cost_function] Best segment = %d -- Volume = %g -- Best position = (%g,%g,%g)\n",\
                                    best->id,\
                                    minimum_volume,\
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
            // Volume
            double volume = calc_tree_volume(the_network);

            // Angle
            struct point *src = iconn->value->src->value;
            struct point *dest = iconn->value->dest->value;

            double middle_pos[3];
            calc_middle_point_segment(iconn,middle_pos);

            double u[3], v[3];
            build_unitary_vector(u,middle_pos[0],middle_pos[1],middle_pos[2],\
                                dest->x,dest->y,dest->z);
            build_unitary_vector(v,middle_pos[0],middle_pos[1],middle_pos[2],\
                                new_pos[0],new_pos[1],new_pos[2]);

            double dist = euclidean_norm(new_pos[0],new_pos[1],new_pos[2],\
                                        middle_pos[0],middle_pos[1],middle_pos[2]);

            double degrees = calc_angle_between_vectors(u,v);
            if (volume < minimum_volume && degrees > degrees_limit)
            {
                minimum_volume = volume;
                best = iconn;

                printf("[cost_function] Best segment = %d -- Volume = %g\n",best->id,minimum_volume);
            }

            restore_state_tree(the_network,iconn);
        }
    }

    return best;
}

SET_COST_FUNCTION (minimize_tree_volume_with_level_penalty)
{
    struct segment_node *best = NULL;
    double minimum_eval = __DBL_MAX__;

    double Q_perf = the_network->Q_perf;
    double p_perf = the_network->p_perf;
    double p_term = the_network->p_term;
    double delta_p = p_perf - p_term;

    // Get a reference to the local optimization function
    bool using_local_optimization = the_network->using_local_optimization;
    set_local_optimization_function_fn *local_optimization_fn = local_opt_config->function;

    for (uint32_t i = 0; i < feasible_segments.size(); i++)
    {
        struct segment_node *iconn = feasible_segments[i];

        // First we set the bifurcation point on the middle of the segment
        struct segment_node *inew = build_segment(the_network,local_opt_config,iconn->id,new_pos);

        struct segment_node *ibiff = iconn->value->parent;

        //  WITH LOCAL OPTIMIZATION
        if (using_local_optimization)
        {
            // Save the original position of the bifurcation
            double ori_pos[3];
            save_original_bifurcation_position(ibiff,ori_pos);

            // Initialize the best position as middle point of the segment
            double best_pos[3];
            initialize_best_position_as_middle_point(best_pos,ori_pos);

            //printf("[cost_function] Original bifurcation position = (%g,%g,%g)\n",\
                                                        best_pos[0],\
                                                        best_pos[1],\
                                                        best_pos[2]);

            // 1) Check the cost function of the first configuration
            double volume = calc_tree_volume(the_network);
            double eval = calc_segment_custom_function_with_level_penalty(volume,iconn);

            if (eval < minimum_eval)
            {
                minimum_eval = eval;
                best = iconn;

                // The best position of the best segment will be stored inside the 
                // 'local_optimization' structure
                local_opt_config->best_pos[0] = best_pos[0];
                local_opt_config->best_pos[1] = best_pos[1];
                local_opt_config->best_pos[2] = best_pos[2];

                printf("[cost_function] Best segment = %d -- Eval = %g -- Best position = (%g,%g,%g)\n",\
                                best->id,\
                                minimum_eval,\
                                local_opt_config->best_pos[0],\
                                local_opt_config->best_pos[1],\
                                local_opt_config->best_pos[2]);

            }

            // 2) Now, call the local optimization function and fill the 'test_positions' array
            std::vector<struct point> test_positions;
            local_optimization_fn(iconn,ibiff,inew,test_positions);

            for (uint32_t j = 0; j < test_positions.size(); j++)
            {
                // Change the position of the bifurcation point 
                double new_pos[3];
                new_pos[0] = test_positions[j].x;
                new_pos[1] = test_positions[j].y;
                new_pos[2] = test_positions[j].z;
                move_bifurcation_location(iconn,ibiff,inew,new_pos);

                rescale_tree(ibiff,iconn,inew,Q_perf,delta_p,the_network->num_terminals);

                recalculate_radius(the_network);

                double volume = calc_tree_volume(the_network);
                double eval = calc_segment_custom_function_with_level_penalty(volume,iconn);

                if (eval < minimum_eval)
                {
                    minimum_eval = eval;
                    best = iconn;

                    // The best position of the best segment will be stored inside the 'local_optimization' structure
                    local_opt_config->best_pos[0] = new_pos[0];
                    local_opt_config->best_pos[1] = new_pos[1];
                    local_opt_config->best_pos[2] = new_pos[2];

                    printf("[cost_function] Best segment = %d -- Eval = %g -- Best position = (%g,%g,%g)\n",\
                                    best->id,\
                                    minimum_eval,\
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
            double volume = calc_tree_volume(the_network);
            double eval = calc_segment_custom_function_with_level_penalty(volume,iconn);

            if (eval < minimum_eval)
            {
                minimum_eval = eval;
                best = iconn;

                printf("[cost_function] Best segment = %d -- Eval = %g\n",best->id,minimum_eval);
            }

            restore_state_tree(the_network,iconn);
        }
    }

    return best;
}
