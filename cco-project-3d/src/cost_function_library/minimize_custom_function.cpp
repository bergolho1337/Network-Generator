#include "minimize_custom_function.h"

SET_COST_FUNCTION (minimize_custom_function)
{
    FILE *log_file = the_network->log_file;

    struct segment_node *best = NULL;
    double minimum_eval = __DBL_MAX__;

    double Q_perf = the_network->Q_perf;
    double p_perf = the_network->p_perf;
    double p_term = the_network->p_term;
    double delta_p = p_perf - p_term;

    // Get a reference to the local optimization function
    bool using_local_optimization = the_network->using_local_optimization;
    set_local_optimization_function_fn *local_optimization_fn = local_opt_config->function;

    double beta;
    get_parameter_value_from_map(config->params,"beta",&beta);

    double alpha;
    get_parameter_value_from_map(config->params,"alpha",&alpha);

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
            save_original_bifurcation_position(ibiff,ori_pos);

            // Initialize the best position as middle point of the segment
            double best_pos[3];
            initialize_best_position_as_middle_point(best_pos,ori_pos);

            //printf("[cost_function] Original bifurcation position = (%g,%g,%g)\n",\
                                                        best_pos[0],\
                                                        best_pos[1],\
                                                        best_pos[2]);

            // 1) Check the cost function of the first configuration
            double eval = calc_custom_function(the_network,beta,alpha);

            // Collision detection: Check if the new segment 'inew' collides with any other segment from the network a part from the 'iconn'
            struct segment_list *s_list = the_network->segment_list;
            bool point_is_not_ok = has_collision(s_list,iconn,ibiff,inew,log_file);

            if (eval < minimum_eval && !point_is_not_ok)
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
                double new_biff_pos[3];
                new_biff_pos[0] = test_positions[j].x;
                new_biff_pos[1] = test_positions[j].y;
                new_biff_pos[2] = test_positions[j].z;
                move_bifurcation_location(iconn,ibiff,inew,new_biff_pos);

                rescale_tree(ibiff,iconn,inew,Q_perf,delta_p,the_network->num_terminals,the_network->using_only_murray_law);

                recalculate_radius(the_network);

                double eval = calc_custom_function(the_network,beta,alpha);

                // Collision detection: Check if the new segment 'inew' collides with any other segment from the network a part from the 'iconn' and 'ibiff'
                s_list = the_network->segment_list;
                point_is_not_ok = has_collision(s_list,iconn,ibiff,inew,log_file);

                if (eval < minimum_eval && !point_is_not_ok)
                {
                    minimum_eval = eval;
                    best = iconn;

                    // The best position of the best segment will be stored inside the 'local_optimization' structure
                    local_opt_config->best_pos[0] = new_biff_pos[0];
                    local_opt_config->best_pos[1] = new_biff_pos[1];
                    local_opt_config->best_pos[2] = new_biff_pos[2];

                    printf("[cost_function] Best segment = %d -- Eval = %g -- Best position = (%g,%g,%g) [%u]\n",\
                                    best->id,\
                                    minimum_eval,\
                                    local_opt_config->best_pos[0],\
                                    local_opt_config->best_pos[1],\
                                    local_opt_config->best_pos[2],j);
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
            double eval = calc_custom_function(the_network,beta,alpha);

            if (eval < minimum_eval)
            {
                minimum_eval = eval;
                best = iconn;

                printf("[cost_function] Best segment = %d -- Eval = %g\n",best->id,minimum_eval);
            }

            restore_state_tree(the_network,iconn);
        }

    }

    if (using_local_optimization)
    {
        // Set off the 'first_call' flag
        local_opt_config->first_call = false;
    }

    return best;
}

/*
SET_COST_FUNCTION (minimize_custom_function_with_angle_restriction)
{
    FILE *log_file = the_network->log_file;

    struct segment_node *best = NULL;
    double minimum_eval = __DBL_MAX__;

    double Q_perf = the_network->Q_perf;
    double p_perf = the_network->p_perf;
    double p_term = the_network->p_term;
    double delta_p = p_perf - p_term;

    // Get a reference to the local optimization function
    bool using_local_optimization = the_network->using_local_optimization;
    set_local_optimization_function_fn *local_optimization_fn = local_opt_config->function;

    double beta;
    get_parameter_value_from_map(config->params,"beta",&beta);
    double alpha;
    get_parameter_value_from_map(config->params,"alpha",&alpha);
    double min_degrees_limit;
    get_parameter_value_from_map(config->params,"min_degrees_limit",&min_degrees_limit);
    double max_degrees_limit;
    get_parameter_value_from_map(config->params,"max_degrees_limit",&max_degrees_limit);

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
            save_original_bifurcation_position(ibiff,ori_pos);

            // Initialize the best position as middle point of the segment
            double best_pos[3];
            initialize_best_position_as_middle_point(best_pos,ori_pos);

            //printf("[cost_function] Original bifurcation position = (%g,%g,%g)\n",\
                                                        best_pos[0],\
                                                        best_pos[1],\
                                                        best_pos[2]);

            // 1) Check the cost function of the first configuration
            double eval = calc_custom_function(the_network,beta,alpha);

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

            // Collision detection: Check if the new segment 'inew' collides with any other segment from the network a part from the 'iconn'
            struct segment_list *s_list = the_network->segment_list;
            bool point_is_not_ok = has_collision(s_list,iconn,ibiff,inew,log_file);

            if (eval < minimum_eval && degrees > min_degrees_limit && degrees < max_degrees_limit && !point_is_not_ok)
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
                double new_biff_pos[3];
                new_biff_pos[0] = test_positions[j].x;
                new_biff_pos[1] = test_positions[j].y;
                new_biff_pos[2] = test_positions[j].z;
                move_bifurcation_location(iconn,ibiff,inew,new_biff_pos);

                rescale_tree(ibiff,iconn,inew,Q_perf,delta_p,the_network->num_terminals,the_network->using_only_murray_law);

                recalculate_radius(the_network);

                double eval = calc_custom_function(the_network,beta,alpha);

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

                // Collision detection: Check if the new segment 'inew' collides with any other segment from the network a part from the 'iconn' and 'ibiff'
                s_list = the_network->segment_list;
                point_is_not_ok = has_collision(s_list,iconn,ibiff,inew,log_file);

                if (eval < minimum_eval && degrees > min_degrees_limit && degrees < max_degrees_limit && !point_is_not_ok)
                {
                    minimum_eval = eval;
                    best = iconn;

                    // The best position of the best segment will be stored inside the 'local_optimization' structure
                    local_opt_config->best_pos[0] = new_biff_pos[0];
                    local_opt_config->best_pos[1] = new_biff_pos[1];
                    local_opt_config->best_pos[2] = new_biff_pos[2];

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
            double eval = calc_custom_function(the_network,beta,alpha);

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

            if (eval < minimum_eval && degrees > min_degrees_limit && degrees < max_degrees_limit)
            {
                minimum_eval = eval;
                best = iconn;

                printf("[cost_function] Best segment = %d -- Eval = %g\n",best->id,minimum_eval);
            }

            restore_state_tree(the_network,iconn);
        }

    }

    if (using_local_optimization)
    {
        // Set off the 'first_call' flag
        local_opt_config->first_call = false;
    }

    return best;
}

// TODO: Need improvement ...
SET_COST_FUNCTION (minimize_custom_function_with_level_penalty)
{
    FILE *log_file = the_network->log_file;

    struct segment_node *best = NULL;
    double minimum_eval = __DBL_MAX__;

    double Q_perf = the_network->Q_perf;
    double p_perf = the_network->p_perf;
    double p_term = the_network->p_term;
    double delta_p = p_perf - p_term;

    // Get a reference to the local optimization function
    bool using_local_optimization = the_network->using_local_optimization;
    set_local_optimization_function_fn *local_optimization_fn = local_opt_config->function;

    double beta;
    get_parameter_value_from_map(config->params,"beta",&beta);

    double alpha;
    get_parameter_value_from_map(config->params,"alpha",&alpha);

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
            save_original_bifurcation_position(ibiff,ori_pos);

            // Initialize the best position as middle point of the segment
            double best_pos[3];
            initialize_best_position_as_middle_point(best_pos,ori_pos);

            //printf("[cost_function] Original bifurcation position = (%g,%g,%g)\n",\
                                                        best_pos[0],\
                                                        best_pos[1],\
                                                        best_pos[2]);

            // 1) Check the cost function of the first configuration
            double eval = calc_custom_function(the_network,beta,alpha);
            eval = calc_segment_custom_function_with_level_penalty(eval,iconn);

            // Collision detection: Check if the new segment 'inew' collides with any other segment from the network a part from the 'iconn'
            struct segment_list *s_list = the_network->segment_list;
            bool point_is_not_ok = has_collision(s_list,iconn,ibiff,inew,log_file);

            if (eval < minimum_eval && !point_is_not_ok)
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
                double new_biff_pos[3];
                new_biff_pos[0] = test_positions[j].x;
                new_biff_pos[1] = test_positions[j].y;
                new_biff_pos[2] = test_positions[j].z;
                move_bifurcation_location(iconn,ibiff,inew,new_biff_pos);

                rescale_tree(ibiff,iconn,inew,Q_perf,delta_p,the_network->num_terminals,the_network->using_only_murray_law);

                recalculate_radius(the_network);

                double eval = calc_custom_function(the_network,beta,alpha);
                eval = calc_segment_custom_function_with_level_penalty(eval,iconn);

                // Collision detection: Check if the new segment 'inew' collides with any other segment from the network a part from the 'iconn' and 'ibiff'
                s_list = the_network->segment_list;
                point_is_not_ok = has_collision(s_list,iconn,ibiff,inew,log_file);

                if (eval < minimum_eval && !point_is_not_ok)
                {
                    minimum_eval = eval;
                    best = iconn;

                    // The best position of the best segment will be stored inside the 'local_optimization' structure
                    local_opt_config->best_pos[0] = new_biff_pos[0];
                    local_opt_config->best_pos[1] = new_biff_pos[1];
                    local_opt_config->best_pos[2] = new_biff_pos[2];

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
            double eval = calc_custom_function(the_network,beta,alpha);
            eval = calc_segment_custom_function_with_level_penalty(eval,iconn);

            if (eval < minimum_eval)
            {
                minimum_eval = eval;
                best = iconn;

                printf("[cost_function] Best segment = %d -- Eval = %g\n",best->id,minimum_eval);
            }

            restore_state_tree(the_network,iconn);
        }
    }

    if (using_local_optimization)
    {
        // Set off the 'first_call' flag
        local_opt_config->first_call = false;
    }

    return best;
}
*/