#include "minimize_custom_function.h"

SET_COST_FUNCTION (minimize_custom_function)
{
    struct segment_node *best = NULL;
    double minimum_eval = __DBL_MAX__;

    double I_in = the_network->I_in;
    double V_in = the_network->V_in;
    double V_out = the_network->V_out;
    double delta_v = V_in - V_out;

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

                rescale_tree(ibiff,iconn,inew,I_in,delta_v,the_network->num_terminals);

                recalculate_radius(the_network);

                double eval = calc_custom_function(the_network,beta,alpha);

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

    return best;
}

// TODO: Need improvement ...
SET_COST_FUNCTION (minimize_custom_function_with_level_penalty)
{
    struct segment_node *best = NULL;
    double minimum_eval = __DBL_MAX__;

    double I_in = the_network->I_in;
    double V_in = the_network->V_in;
    double V_out = the_network->V_out;
    double delta_v = V_in - V_out;

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

                rescale_tree(ibiff,iconn,inew,I_in,delta_v,the_network->num_terminals);

                recalculate_radius(the_network);

                double eval = calc_custom_function(the_network,beta,alpha);
                eval = calc_segment_custom_function_with_level_penalty(eval,iconn);

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

    return best;
}