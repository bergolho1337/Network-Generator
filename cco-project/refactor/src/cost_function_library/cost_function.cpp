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
    double closest_dist = DBL_MAX;

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