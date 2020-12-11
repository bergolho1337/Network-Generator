#include "cost_function.h"

#include "../cco/cco.h"

void CustomFunction::init_parameters (CostFunctionConfig *cost_function_config)
{
    beta = 1.0;
    cost_function_config->get_parameter_value_from_map("beta",&beta);
    alpha = 0.0;
    cost_function_config->get_parameter_value_from_map("alpha",&alpha);
    min_degrees_limit = 1.0;
    cost_function_config->get_parameter_value_from_map("min_degrees_limit",&min_degrees_limit);
    max_degrees_limit = 180.0;
    cost_function_config->get_parameter_value_from_map("max_degrees_limit",&max_degrees_limit);
    min_segment_length = 0.0001;
    cost_function_config->get_parameter_value_from_map("min_segment_length",&min_segment_length);
    max_segment_length = 0.01;
    cost_function_config->get_parameter_value_from_map("max_segment_length",&max_segment_length);
}

bool CustomFunction::check_restrictions (CCO_Network *the_network, Segment *iconn, Segment *ibiff, Segment *inew)
{
    // Check bifurcation size and angle
    bool has_angle_requirement = check_angle_restriction(iconn,inew);
    
    // Check segment sizes
    bool has_minimum_segment_size = check_minimum_segment_size(iconn,ibiff,inew);
    bool has_maximum_segment_size = check_maximum_segment_size(iconn,ibiff,inew);

    // Collision detection: Check if the new segment 'inew' collides with any other segment from the network a part from the 'iconn'
    bool has_segment_segment_collision = check_collision(the_network,iconn,ibiff,inew);
    //bool has_segment_triangle_collision = has_intersect_obstacle(inew,obstacle_faces); 

    return has_segment_segment_collision ||\
           !has_minimum_segment_size ||\
           !has_maximum_segment_size ||\
           !has_angle_requirement;
}

Segment* CustomFunction::eval (CCO_Network *the_network,\
                        CostFunctionConfig *cost_function_config,\
                        LocalOptimizationConfig *local_opt_config,\
                        std::vector<Segment*> feasible_segments,\
                        Point *new_term)
{
    FILE *log_file = the_network->log_file;

    Segment *best = NULL;
    double minimum_eval = __DBL_MAX__;

    bool using_local_optimization = the_network->using_local_optimization;
    LocalOptimization *local_opt_fn = the_network->local_opt_fn;

    // Get cost function parameters or use the default ones
    init_parameters(cost_function_config);
    
    double u[3], v[3];
    double angle_degrees;

    for (uint32_t i = 0; i < feasible_segments.size(); i++)
    {
        Segment *iconn = feasible_segments[i];
        Segment *inew = the_network->build_segment(local_opt_config,iconn->id,new_term);
        Segment *ibiff = iconn->parent;

        if (the_network->using_local_optimization)
        {
            // Save the original position of the bifurcation
            double ori_pos[3];
            local_opt_config->save_original_bifurcation_position(ibiff,ori_pos);

            // Initialize the best position as middle point of the segment
            double best_pos[3];
            local_opt_config->initialize_best_position_as_middle_point(best_pos,ori_pos);

            // [EVALUATE COST FUNCTION]
            double eval = calc_custom_function(the_network);

            // [RESTRICTION SECTION]
            bool point_is_not_ok = check_restrictions(the_network,iconn,ibiff,inew);

            if (eval < minimum_eval && !point_is_not_ok)
            {
                minimum_eval = eval;
                best = iconn;

                // The best position of the best segment will be stored inside the 
                // 'local_optimization' structure
                memcpy(local_opt_config->best_pos,best_pos,sizeof(double)*3);

                printf("[cost_function] Best segment = %d -- Eval = %g -- Best position = (%g,%g,%g)\n",\
                                best->id,\
                                minimum_eval,\
                                local_opt_config->best_pos[0],\
                                local_opt_config->best_pos[1],\
                                local_opt_config->best_pos[2]);
            }

            // 2) Now, call the local optimization function and fill the 'test_positions' array
            std::vector<Point*> test_positions;
            local_opt_fn->optimize(iconn,ibiff,inew,test_positions);

            for (uint32_t j = 0; j < test_positions.size(); j++)
            {
                // Change the position of the bifurcation point 
                local_opt_fn->move_bifurcation_location(iconn,ibiff,inew,test_positions[j]);

                // Rescale the tree using the current configuration
                the_network->rescale_tree(ibiff,iconn,inew);
                the_network->recalculate_radius();
                the_network->recalculate_length();

                // [EVALUATE COST FUNCTION]
                double eval = calc_custom_function(the_network);

                // [RESTRICTION SECTION]
                bool point_is_not_ok = check_restrictions(the_network,iconn,ibiff,inew);

                if (eval < minimum_eval && !point_is_not_ok)
                {
                    minimum_eval = eval;
                    best = iconn;

                    // The best position of the best segment will be stored inside the 'local_optimization' structure
                    local_opt_config->update_bifurcation_position(test_positions[j]);

                    printf("[cost_function] Best segment = %d -- Eval = %g -- Best position = (%g,%g,%g)\n",\
                                    best->id,\
                                    minimum_eval,\
                                    local_opt_config->best_pos[0],\
                                    local_opt_config->best_pos[1],\
                                    local_opt_config->best_pos[2]);
                }
            }

            // Move the bifurcation to the original position
            local_opt_fn->move_bifurcation_location(iconn,ibiff,inew,ori_pos);

            // Restore the tree configuration before the 
            the_network->restore_state_tree(iconn);

            // Clear the array for the next segment
            local_opt_fn->free_test_positions(test_positions);
            test_positions.clear();
        }
        // NO LOCAL OPTIMIZATION
        else
        {
            // [EVALUATE COST FUNCTION]
            double eval = calc_custom_function(the_network);

            // [RESTRICTION SECTION]
            bool point_is_not_ok = check_restrictions(the_network,iconn,ibiff,inew);

            if (eval < minimum_eval && !point_is_not_ok)
            {
                minimum_eval = eval;
                best = iconn;

                printf("[cost_function] Best segment = %d -- Eval = %g\n",best->id,minimum_eval);
            }

            the_network->restore_state_tree(iconn);
        }
    }
    
    if (the_network->using_local_optimization)
    {
        // Set off the 'first_call' flag
        local_opt_config->first_call = false;
    }

    return best;
}

double CustomFunction::calc_custom_function (CCO_Network *the_network)
{
    double result = 0.0;
    for (uint32_t i = 0; i < the_network->segment_list.size(); i++)
    {
        Segment *cur_segment = the_network->segment_list[i];

        result += calc_segment_custom_function(cur_segment);
    }
    return result;
}

double CustomFunction::calc_segment_custom_function (Segment *s)
{
    double l = s->length;
    double r = s->radius;

    double length = euclidean_norm(s->src->x,s->src->y,s->src->z,\
                                s->dest->x,s->dest->y,s->dest->z);
    
    return pow(l,beta) * pow(r,alpha);
}

bool CustomFunction::check_angle_restriction (Segment *iconn, Segment *inew)
{
    double u[3], v[3];
    
    iconn->calc_unitary_vector(u);
    inew->calc_unitary_vector(v);
    double angle = calc_angle_between_vectors(u,v);

    return (angle > min_degrees_limit && angle < max_degrees_limit);
}

bool CustomFunction::check_minimum_segment_size (Segment *iconn, Segment *ibiff, Segment *inew)
{
    return (iconn->length > min_segment_length && ibiff->length > min_segment_length && inew->length > min_segment_length);
}

bool CustomFunction::check_maximum_segment_size (Segment *iconn, Segment *ibiff, Segment *inew)
{
    return (iconn->length < max_segment_length && ibiff->length < max_segment_length && inew->length < max_segment_length);
}

bool CustomFunction::check_collision (CCO_Network *the_network, Segment *iconn, Segment *ibiff, Segment *inew)
{
    // Get the reference to the points from the 'inew' segment
    Point *src_inew = inew->src;
    Point *dest_inew = inew->dest;

    // Get the reference to the points from the 'ibiff' segment
    Point *src_ibiff = ibiff->src;
    Point *dest_ibiff = ibiff->dest;

    for (uint32_t i = 0; i < the_network->segment_list.size(); i++)
    {
        Segment *cur_segment = the_network->segment_list[i];

        // We avoid both the 'iconn' and 'inew' segments from the search
        if (cur_segment->id != iconn->id && cur_segment->id != inew->id && cur_segment->id != ibiff->id)
        {
            // Get the reference to the points of the current segment
            Point *src = cur_segment->src;
            Point *dest = cur_segment->dest;

            bool intersect = collision_detection(src_inew->x,src_inew->y,src_inew->z,\
                                            dest_inew->x,dest_inew->y,dest_inew->z,\
                                            src->x,src->y,src->z,\
                                            dest->x,dest->y,dest->z);

            if (intersect)
            {
                //printf("\t[-] ERROR! Intersection between segments!\n");
                fprintf(the_network->log_file,"\t[-] ERROR! Intersection with segment %d !\n",cur_segment->id);
                return true;
            }
        }
    }
    return false;
}




void ActivationTimeFunction::init_parameters (CostFunctionConfig *cost_function_config)
{
    // Parameters of the function goes here ...
}

bool ActivationTimeFunction::check_restrictions (CCO_Network *the_network, Segment *iconn, Segment *ibiff, Segment *inew)
{
    return true;    
}

Segment* ActivationTimeFunction::eval (CCO_Network *the_network,\
                        CostFunctionConfig *cost_function_config,\
                        LocalOptimizationConfig *local_opt_config,\
                        std::vector<Segment*> feasible_segments,\
                        Point *new_term)
{
    printf("Activation time function\n");
    return NULL;
}