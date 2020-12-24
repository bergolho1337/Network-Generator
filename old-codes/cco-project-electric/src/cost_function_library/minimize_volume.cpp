#include "minimize_volume.h"

// OK
SET_COST_FUNCTION (minimize_tree_volume_default)
{
    FILE *log_file = the_network->log_file;

    struct segment_node *best = NULL;
    double minimum_volume = __DBL_MAX__;

    double I_in = the_network->I_in;
    double V_in = the_network->V_in;
    double V_out = the_network->V_out;
    double delta_v = V_in - V_out;

    // Get a reference to the local optimization function
    set_local_optimization_function_fn *local_optimization_fn = local_opt_config->function;
    bool using_local_optimization = the_network->using_local_optimization;

    // For each feasible segment try to make a connection
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

            // Collision detection: Check if the new segment 'inew' collides with any other segment from the network a part from the 'iconn'
            struct segment_list *s_list = the_network->segment_list;
            bool point_is_not_ok = has_collision(s_list,iconn,ibiff,inew,log_file);

            if (volume < minimum_volume && !point_is_not_ok)
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
                double new_biff_pos[3];
                new_biff_pos[0] = test_positions[j].x;
                new_biff_pos[1] = test_positions[j].y;
                new_biff_pos[2] = test_positions[j].z;
                move_bifurcation_location(iconn,ibiff,inew,new_biff_pos);

                rescale_tree(ibiff,iconn,inew,I_in,delta_v,the_network->num_terminals);

                recalculate_radius(the_network);

                double volume = calc_tree_volume(the_network);

                // Collision detection: Check if the new segment 'inew' collides with any other segment from the network a part from the 'iconn' and 'ibiff'
                s_list = the_network->segment_list;
                point_is_not_ok = has_collision(s_list,iconn,ibiff,inew,log_file);

                if (volume < minimum_volume && !point_is_not_ok)
                {
                    minimum_volume = volume;
                    best = iconn;

                    // The best position of the best segment will be stored inside the 'local_optimization' structure
                    local_opt_config->best_pos[0] = new_biff_pos[0];
                    local_opt_config->best_pos[1] = new_biff_pos[1];
                    local_opt_config->best_pos[2] = new_biff_pos[2];

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

    if (using_local_optimization)
    {
        // Set off the 'first_call' flag
        local_opt_config->first_call = false;
    }

    return best;
}

// OK
SET_COST_FUNCTION (minimize_tree_volume_with_angle_restriction)
{
    FILE *log_file = the_network->log_file;

    struct segment_node *best = NULL;
    double minimum_volume = __DBL_MAX__;

    double I_in = the_network->I_in;
    double V_in = the_network->V_in;
    double V_out = the_network->V_out;
    double delta_v = V_in - V_out;

    // Get a reference to the local optimization function
    bool using_local_optimization = the_network->using_local_optimization;
    set_local_optimization_function_fn *local_optimization_fn = local_opt_config->function;

    // Get the degree limit from the user input
    double min_degrees_limit;
    get_parameter_value_from_map(config->params,"min_degrees_limit",&min_degrees_limit);
    double max_degrees_limit;
    get_parameter_value_from_map(config->params,"max_degrees_limit",&max_degrees_limit);

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

            // Collision detection: Check if the new segment 'inew' collides with any other segment from the network a part from the 'iconn'
            struct segment_list *s_list = the_network->segment_list;
            bool point_is_not_ok = has_collision(s_list,iconn,ibiff,inew,log_file);

            if (volume < minimum_volume && degrees > min_degrees_limit && degrees < max_degrees_limit && !point_is_not_ok)
            {
                minimum_volume = volume;
                best = iconn;

                // The best position of the best segment will be stored inside the 
                // 'local_optimization' structure
                local_opt_config->best_pos[0] = best_pos[0];
                local_opt_config->best_pos[1] = best_pos[1];
                local_opt_config->best_pos[2] = best_pos[2];

                printf("[cost_function] Best segment = %d -- Volume = %g -- Best position = (%g,%g,%g) -- Degree = %g\n",\
                                best->id,\
                                minimum_volume,\
                                local_opt_config->best_pos[0],\
                                local_opt_config->best_pos[1],\
                                local_opt_config->best_pos[2],\
                                degrees);

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

                rescale_tree(ibiff,iconn,inew,I_in,delta_v,the_network->num_terminals);

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

                // Collision detection: Check if the new segment 'inew' collides with any other segment from the network a part from the 'iconn' and 'ibiff'
                s_list = the_network->segment_list;
                point_is_not_ok = has_collision(s_list,iconn,ibiff,inew,log_file);

                if (volume < minimum_volume && degrees > min_degrees_limit && degrees < max_degrees_limit && !point_is_not_ok)
                {
                    minimum_volume = volume;
                    best = iconn;

                    // The best position of the best segment will be stored inside the 'local_optimization' structure
                    local_opt_config->best_pos[0] = new_biff_pos[0];
                    local_opt_config->best_pos[1] = new_biff_pos[1];
                    local_opt_config->best_pos[2] = new_biff_pos[2];

                    printf("[cost_function] Best segment = %d -- Volume = %g -- Best position = (%g,%g,%g) -- Degree = %g\n",\
                                    best->id,\
                                    minimum_volume,\
                                    local_opt_config->best_pos[0],\
                                    local_opt_config->best_pos[1],\
                                    local_opt_config->best_pos[2],\
                                    degrees);
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
            if (volume < minimum_volume && degrees > min_degrees_limit && degrees < max_degrees_limit)
            {
                minimum_volume = volume;
                best = iconn;

                printf("[cost_function] Best segment = %d -- Volume = %g\n",best->id,minimum_volume);
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

// OK
SET_COST_FUNCTION (minimize_tree_volume_with_length_restriction)
{
    FILE *log_file = the_network->log_file;

    struct segment_node *best = NULL;
    double minimum_volume = __DBL_MAX__;

    double I_in = the_network->I_in;
    double V_in = the_network->V_in;
    double V_out = the_network->V_out;
    double delta_v = V_in - V_out;

    // Get a reference to the local optimization function
    bool using_local_optimization = the_network->using_local_optimization;
    set_local_optimization_function_fn *local_optimization_fn = local_opt_config->function;

    // Get the length limit from the user input
    double length_limit;
    get_parameter_value_from_map(config->params,"length_limit",&length_limit);

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

            // 1.2) Calculate new segment length
            struct point *src = iconn->value->src->value;
            struct point *dest = iconn->value->dest->value;

            double middle_pos[3];
            calc_middle_point_segment(iconn,middle_pos);

            double length = euclidean_norm(new_pos[0],new_pos[1],new_pos[2],\
                                        middle_pos[0],middle_pos[1],middle_pos[2]);

            // Collision detection: Check if the new segment 'inew' collides with any other segment from the network a part from the 'iconn'
            struct segment_list *s_list = the_network->segment_list;
            bool point_is_not_ok = has_collision(s_list,iconn,ibiff,inew,log_file);

            if (volume < minimum_volume && length > length_limit && !point_is_not_ok)
            {
                minimum_volume = volume;
                best = iconn;

                // The best position of the best segment will be stored inside the 
                // 'local_optimization' structure
                local_opt_config->best_pos[0] = best_pos[0];
                local_opt_config->best_pos[1] = best_pos[1];
                local_opt_config->best_pos[2] = best_pos[2];

                printf("[cost_function] Best segment = %d -- Volume = %g -- Length = %g -- Best position = (%g,%g,%g)\n",\
                                best->id,\
                                minimum_volume,\
                                length,\
                                local_opt_config->best_pos[0],\
                                local_opt_config->best_pos[1],\
                                local_opt_config->best_pos[2]);

            }

            // 2) Now, call the local optimization function and fill the 'test_positions' array
            std::vector<struct point> test_positions;
            local_optimization_fn(iconn,ibiff,inew,test_positions);

            for (uint32_t j = 0; j < test_positions.size(); j++)
            {
                // 2.1) Change the position of the bifurcation point 
                double new_biff_pos[3];
                new_biff_pos[0] = test_positions[j].x;
                new_biff_pos[1] = test_positions[j].y;
                new_biff_pos[2] = test_positions[j].z;
                move_bifurcation_location(iconn,ibiff,inew,new_biff_pos);

                rescale_tree(ibiff,iconn,inew,I_in,delta_v,the_network->num_terminals);

                recalculate_radius(the_network);

                // 2.2) Calculate the volume of the tree on the new configuration
                double volume = calc_tree_volume(the_network);

                // 2.3) Calculate the length of the new segment
                double middle_pos[3];
                calc_middle_point_segment(iconn,middle_pos);

                double length = euclidean_norm(new_pos[0],new_pos[1],new_pos[2],\
                                            middle_pos[0],middle_pos[1],middle_pos[2]);

                // Collision detection: Check if the new segment 'inew' collides with any other segment from the network a part from the 'iconn' and 'ibiff'
                s_list = the_network->segment_list;
                point_is_not_ok = has_collision(s_list,iconn,ibiff,inew,log_file);

                if (volume < minimum_volume && length > length_limit && !point_is_not_ok)
                {
                    minimum_volume = volume;
                    best = iconn;

                    // The best position of the best segment will be stored inside the 'local_optimization' structure
                    local_opt_config->best_pos[0] = new_biff_pos[0];
                    local_opt_config->best_pos[1] = new_biff_pos[1];
                    local_opt_config->best_pos[2] = new_biff_pos[2];

                    printf("[cost_function] Best segment = %d -- Volume = %g -- Length = %g -- Best position = (%g,%g,%g)\n",\
                                    best->id,\
                                    minimum_volume,\
                                    length,\
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
            // Calculate the volume of the tree
            double volume = calc_tree_volume(the_network);

            // Calculate the length of the new segment
            double middle_pos[3];
            calc_middle_point_segment(iconn,middle_pos);

            double length = euclidean_norm(new_pos[0],new_pos[1],new_pos[2],\
                                        middle_pos[0],middle_pos[1],middle_pos[2]);

            if (volume < minimum_volume && length > length_limit)
            {
                minimum_volume = volume;
                best = iconn;

                printf("[cost_function] Best segment = %d -- Length = %g -- Volume = %g\n",\
                                    best->id,\
                                    length,\
                                    minimum_volume);
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

// OK
SET_COST_FUNCTION (minimize_tree_volume_with_angle_and_length_restriction)
{
    FILE *log_file = the_network->log_file;

    struct segment_node *best = NULL;
    double minimum_volume = __DBL_MAX__;

    double I_in = the_network->I_in;
    double V_in = the_network->V_in;
    double V_out = the_network->V_out;
    double delta_v = V_in - V_out;

    // Get a reference to the local optimization function
    bool using_local_optimization = the_network->using_local_optimization;
    set_local_optimization_function_fn *local_optimization_fn = local_opt_config->function;

    // Get the degree and length limit from the user input
    double min_degrees_limit;
    get_parameter_value_from_map(config->params,"min_degrees_limit",&min_degrees_limit);
    double max_degrees_limit;
    get_parameter_value_from_map(config->params,"max_degrees_limit",&max_degrees_limit);
    double length_limit;
    get_parameter_value_from_map(config->params,"length_limit",&length_limit);

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

            // 3) Calculate new segment length
            double length = dist;

            // Collision detection: Check if the new segment 'inew' collides with any other segment from the network a part from the 'iconn'
            struct segment_list *s_list = the_network->segment_list;
            bool point_is_not_ok = has_collision(s_list,iconn,ibiff,inew,log_file);

            if (volume < minimum_volume && degrees > min_degrees_limit && degrees < max_degrees_limit && length > length_limit && !point_is_not_ok)
            {
                minimum_volume = volume;
                best = iconn;

                // The best position of the best segment will be stored inside the 
                // 'local_optimization' structure
                local_opt_config->best_pos[0] = best_pos[0];
                local_opt_config->best_pos[1] = best_pos[1];
                local_opt_config->best_pos[2] = best_pos[2];

                printf("[cost_function] Best segment = %d -- Volume = %g -- Degrees = %g -- Length = %g -- Best position = (%g,%g,%g)\n",\
                                best->id,\
                                minimum_volume,\
                                degrees,\
                                length,\
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

                rescale_tree(ibiff,iconn,inew,I_in,delta_v,the_network->num_terminals);

                recalculate_radius(the_network);

                // 2.1) Calculate tree volume
                double volume = calc_tree_volume(the_network);

                // 2.2) Calculate bifurcation angle
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

                // 2.3) Calculate new segment length
                double length = dist;

                // Collision detection: Check if the new segment 'inew' collides with any other segment from the network a part from the 'iconn' and 'ibiff'
                s_list = the_network->segment_list;
                point_is_not_ok = has_collision(s_list,iconn,ibiff,inew,log_file);

                if (volume < minimum_volume && degrees > min_degrees_limit && degrees < max_degrees_limit && length > length_limit && !point_is_not_ok)
                {
                    minimum_volume = volume;
                    best = iconn;

                    // The best position of the best segment will be stored inside the 'local_optimization' structure
                    local_opt_config->best_pos[0] = new_biff_pos[0];
                    local_opt_config->best_pos[1] = new_biff_pos[1];
                    local_opt_config->best_pos[2] = new_biff_pos[2];

                    printf("[cost_function] Best segment = %d -- Volume = %g -- Degrees = %g -- Length = %g -- Best position = (%g,%g,%g)\n",\
                                best->id,\
                                minimum_volume,\
                                degrees,\
                                length,\
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

            // Length
            double length = dist;

            if (volume < minimum_volume && degrees > min_degrees_limit && degrees < max_degrees_limit && length > length_limit)
            {
                minimum_volume = volume;
                best = iconn;

                printf("[cost_function] Best segment = %d -- Volume = %g -- Degree = %g -- Length = %g\n",best->id,minimum_volume,degrees,length);
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

// TODO: Need improvements ...
SET_COST_FUNCTION (minimize_tree_volume_with_level_penalty)
{
    FILE *log_file = the_network->log_file;

    struct segment_node *best = NULL;
    double minimum_eval = __DBL_MAX__;

    double I_in = the_network->I_in;
    double V_in = the_network->V_in;
    double V_out = the_network->V_out;
    double delta_v = V_in - V_out;

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

                printf("[cost_function] Best segment = %d -- Eval = %g -- Best position = (%g,%g,%g) -- middle point\n",\
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

                rescale_tree(ibiff,iconn,inew,I_in,delta_v,the_network->num_terminals);

                recalculate_radius(the_network);

                double volume = calc_tree_volume(the_network);
                double eval = calc_segment_custom_function_with_level_penalty(volume,iconn);

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
            double volume = calc_tree_volume(the_network);
            double eval = calc_segment_custom_function_with_level_penalty(volume,iconn);

            // Collision detection
            //struct segment_list *s_list = the_network->segment_list;
            //bool point_is_not_ok = has_collision(s_list,iconn,new_pos,log_file);

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

SET_COST_FUNCTION (minimize_tree_volume_with_level_penalty_and_angle_restriction)
{
    FILE *log_file = the_network->log_file;

    struct segment_node *best = NULL;
    double minimum_eval = __DBL_MAX__;

    double I_in = the_network->I_in;
    double V_in = the_network->V_in;
    double V_out = the_network->V_out;
    double delta_v = V_in - V_out;

    double min_degrees_limit;
    get_parameter_value_from_map(config->params,"min_degrees_limit",&min_degrees_limit);
    double max_degrees_limit;
    get_parameter_value_from_map(config->params,"max_degrees_limit",&max_degrees_limit);

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

                printf("[cost_function] Best segment = %d -- Eval = %g -- Best position = (%g,%g,%g) -- middle point\n",\
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

                rescale_tree(ibiff,iconn,inew,I_in,delta_v,the_network->num_terminals);

                recalculate_radius(the_network);

                // 1) Evaluate the cost function
                double volume = calc_tree_volume(the_network);
                double eval = calc_segment_custom_function_with_level_penalty(volume,iconn);

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

                    printf("[cost_function] Best segment = %d -- Eval = %g -- Degrees = %g -- Best position = (%g,%g,%g)\n",\
                                    best->id,\
                                    minimum_eval,\
                                    degrees,\
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
            // 1) Evaluate the cost function
            double volume = calc_tree_volume(the_network);
            double eval = calc_segment_custom_function_with_level_penalty(volume,iconn);

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

            // Collision detection
            //struct segment_list *s_list = the_network->segment_list;
            //bool point_is_not_ok = has_collision(s_list,iconn,new_pos,log_file);

            if (eval < minimum_eval && degrees > min_degrees_limit && degrees < max_degrees_limit)
            {
                minimum_eval = eval;
                best = iconn;

                printf("[cost_function] Best segment = %d -- Eval = %g -- Degrees = %g\n",best->id,minimum_eval,degrees);
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