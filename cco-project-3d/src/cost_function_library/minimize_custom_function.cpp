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

    // Get cost function parameters
    double beta;
    get_parameter_value_from_map(config->params,"beta",&beta);
    double alpha;
    get_parameter_value_from_map(config->params,"alpha",&alpha);

    // For each feasible segment try to make a connection
    for (uint32_t i = 0; i < feasible_segments.size(); i++)
    {
        struct segment_node *iconn = feasible_segments[i];

        struct segment_node *inew = build_segment(the_network,local_opt_config,iconn->id,new_pos);

        struct segment_node *ibiff = iconn->value->parent;

        // WITH LOCAL OPTIMIZATION
        if (using_local_optimization)
        {
            // Save the original position of the bifurcation
            double ori_pos[3];
            save_original_bifurcation_position(ibiff,ori_pos);

            // Initialize the best position as middle point of the segment
            double best_pos[3];
            initialize_best_position_as_middle_point(best_pos,ori_pos);

            // [EVALUATE COST FUNCTION]
            double eval = calc_custom_function(the_network,beta,alpha);

            // [RESTRICTION SECTION]
            // Bifurcation size
            double iconn_size = calc_segment_size(iconn);
            double ibiff_size = calc_segment_size(ibiff);
            double inew_size = calc_segment_size(inew);
            bool has_minimum_segment_size = has_valid_segment_sizes(iconn_size,ibiff_size,inew_size);

            // Collision detection: Check if the new segment 'inew' collides with any other segment from the network a part from the 'iconn'
            struct segment_list *s_list = the_network->segment_list;
            bool has_segment_segment_collision = has_collision(s_list,iconn,ibiff,inew,log_file);
            bool has_segment_triangle_collision = has_intersect_obstacle(inew,obstacle_faces);

            bool point_is_not_ok = has_segment_segment_collision || has_segment_triangle_collision || !has_minimum_segment_size;

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

            // Now, call the local optimization function and fill the 'test_positions' array
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

                // [EVALUATE COST FUNCTION]
                double eval = calc_custom_function(the_network,beta,alpha);

                // Collision detection: Check if the new segment 'inew' collides with any other segment from the network a part from the 'iconn' and 'ibiff'
                s_list = the_network->segment_list;
                has_segment_segment_collision = has_collision(s_list,iconn,ibiff,inew,log_file);
                has_segment_triangle_collision = has_intersect_obstacle(inew,obstacle_faces);

                // [RESTRICTION ZONE]
                double iconn_size = calc_segment_size(iconn);
                double ibiff_size = calc_segment_size(ibiff);
                double inew_size = calc_segment_size(inew);
                bool has_minimum_segment_size = has_valid_segment_sizes_2(iconn_size,ibiff_size,inew_size);

                point_is_not_ok = has_segment_segment_collision || has_segment_triangle_collision || !has_minimum_segment_size;

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

    // Get cost function parameters
    double beta;
    get_parameter_value_from_map(config->params,"beta",&beta);
    double alpha;
    get_parameter_value_from_map(config->params,"alpha",&alpha);
    double min_degrees_limit;
    get_parameter_value_from_map(config->params,"min_degrees_limit",&min_degrees_limit);
    double max_degrees_limit;
    get_parameter_value_from_map(config->params,"max_degrees_limit",&max_degrees_limit);
    double u[3], v[3];
    double angle_degrees;

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

            // [EVALUATE COST FUNCTION]
            double eval = calc_custom_function(the_network,beta,alpha);

            // [RESTRICTION SECTION]
            // Check bifurcation size and angle
            calc_unitary_vector(iconn,u);
            calc_unitary_vector(inew,v);
            angle_degrees = calc_angle_between_vectors(u,v);
            bool has_angle_requirement = check_angle_restriction(angle_degrees,min_degrees_limit,max_degrees_limit);

            double iconn_size = calc_segment_size(iconn);
            double ibiff_size = calc_segment_size(ibiff);
            double inew_size = calc_segment_size(inew);
            bool has_minimum_segment_size = has_valid_segment_sizes(iconn_size,ibiff_size,inew_size);

            // Collision detection: Check if the new segment 'inew' collides with any other segment from the network a part from the 'iconn'
            struct segment_list *s_list = the_network->segment_list;
            bool has_segment_segment_collision = has_collision(s_list,iconn,ibiff,inew,log_file);
            bool has_segment_triangle_collision = has_intersect_obstacle(inew,obstacle_faces);

            bool point_is_not_ok = has_segment_segment_collision || has_segment_triangle_collision || !has_minimum_segment_size || !has_angle_requirement;
            
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

                // [EVALUATE COST FUNCTION]
                double eval = calc_custom_function(the_network,beta,alpha);

                // [RESTRICTION SECTION]
                // Check bifurcation angle
                calc_unitary_vector(iconn,u);
                calc_unitary_vector(inew,v);
                angle_degrees = calc_angle_between_vectors(u,v);
                has_angle_requirement = check_angle_restriction(angle_degrees,min_degrees_limit,max_degrees_limit);

                double iconn_size = calc_segment_size(iconn);
                double ibiff_size = calc_segment_size(ibiff);
                double inew_size = calc_segment_size(inew);
                has_minimum_segment_size = has_valid_segment_sizes_2(iconn_size,ibiff_size,inew_size);

                // Collision detection: Check if the new segment 'inew' collides with any other segment from the network a part from the 'iconn' and 'ibiff'
                s_list = the_network->segment_list;
                has_segment_segment_collision = has_collision(s_list,iconn,ibiff,inew,log_file);
                has_segment_triangle_collision = has_intersect_obstacle(inew,obstacle_faces);
                
                point_is_not_ok = has_segment_segment_collision || has_segment_triangle_collision || !has_angle_requirement || !has_minimum_segment_size;

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

            calc_unitary_vector(iconn,u);
            calc_unitary_vector(inew,v);
            angle_degrees = calc_angle_between_vectors(u,v);
            bool has_angle_requirement = check_angle_restriction(angle_degrees,min_degrees_limit,max_degrees_limit);

            struct segment_list *s_list = the_network->segment_list;
            bool has_segment_segment_collision = has_collision(s_list,iconn,ibiff,inew,log_file);
            bool has_segment_triangle_collision = has_intersect_obstacle(inew,obstacle_faces);

            bool point_is_not_ok = has_segment_segment_collision || has_segment_triangle_collision || !has_angle_requirement;

            if (eval < minimum_eval && !point_is_not_ok)
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

