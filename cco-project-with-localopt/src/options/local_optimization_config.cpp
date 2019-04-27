#include "local_optimization_config.h"

struct local_optimization_config* new_local_optimization_config ()
{
    struct local_optimization_config *result = (struct local_optimization_config*)malloc(sizeof(struct local_optimization_config));

    result->name = NULL;
    result->handle = NULL;
    
    //result->params = new std::map<std::string,double>();

    return result;
}

void free_local_optimization_config (struct local_optimization_config *config)
{
    //delete config->params;
    if (config->name)
        free(config->name);

    if (config->handle)
        dlclose(config->handle);
    
    free(config);
}

void set_local_optimization_function (struct local_optimization_config *config)
{
    char library_path[MAX_FILENAME_SIZE] = "./shared-libs/libdefault_local_optimization.so";

    config->handle = dlopen(library_path,RTLD_LAZY);
    if (!config->handle) 
    {
        fprintf(stderr,"%s\n",dlerror());
        exit(EXIT_FAILURE);
    }
    else
    {
        fprintf(stdout,"\n[local_optimization] Local optimization library \"%s\" opened with sucess\n",library_path);
    }
    
    char *local_optimiztion_function_name = config->name;

    config->function = (set_local_optimization_function_fn*)dlsym(config->handle,local_optimiztion_function_name);
    if (dlerror() != NULL)  
    {
        fprintf(stderr, "[local_optimization] \"%s\" function not found in the provided local optimization library\n",local_optimiztion_function_name);
        exit(EXIT_FAILURE);
    }
    else
    {
        fprintf(stdout, "[local_optimization] Using the function \"%s\" as local optimization function\n",local_optimiztion_function_name);
    }
}

void print_local_optimization_function_config (struct local_optimization_config *config)
{
    printf("Local optimization function name = \"%s\"\n",config->name);

    /*
    if (!config->params->empty())
    {
        printf("Parameters:\n");
        for (auto it = config->params->begin(); it != config->params->end(); ++it)
            printf("\t%s = %lf\n",it->first.c_str(),it->second);
    }
    */
}

void local_optimization (struct segment_node *iconn, struct segment_node *ibiff, struct segment_node *inew,\
                                std::vector<struct point*> &test_positions)
{
    static const uint32_t NE = 16;
    double delta_e = 1.0 / NE;

    // Get a reference to the 3 vertices surrounding the triangle area
    struct point *g1 = ibiff->value->src->value;
    struct point *g2 = iconn->value->dest->value;
    struct point *g3 = inew->value->dest->value;

    //printf("G1 = (%g,%g,%g)\n",g1->x,g1->y,g1->z);
    //printf("G2 = (%g,%g,%g)\n",g2->x,g2->y,g2->z);
    //printf("G3 = (%g,%g,%g)\n",g3->x,g3->y,g3->z);

    // Build the test points for the local optimization based on the triangle composde by
    // 'iconn', 'ibiff' and 'inew'
    for (uint32_t i = 0; i <= NE; i++)
    {
        for (uint32_t j = 0; j <= NE-i; j++)
        {
            double epsilon = i*delta_e;
            double eta = j*delta_e;

            if (!is_corner(i,j,NE))
            {
                // Build the phi array
                double phi[3];
                phi[0] = 1.0 - epsilon - eta;
                phi[1] = epsilon;
                phi[2] = eta;

                double pos[3];
                pos[0] = (phi[0]*g1->x) + (phi[1]*g2->x) + (phi[2]*g3->x);
                pos[1] = (phi[0]*g1->y) + (phi[1]*g2->y) + (phi[2]*g3->y);
                pos[2] = (phi[0]*g1->z) + (phi[1]*g2->z) + (phi[2]*g3->z);

                struct point *p = new_point(pos);
                test_positions.push_back(p);

            }
        }
    }

}

bool is_corner (const uint32_t i, const uint32_t j, const uint32_t NE)
{
    if ((i == 0 && j == 0) || (i == 0 && j == NE) || (i == NE && j == 0))
        return true;
    else
        return false;
}

void move_bifurcation_location (struct segment_node *iconn, struct segment_node *ibiff, struct segment_node *inew,\
                            const double pos[])
{
    ibiff->value->dest->value->x = pos[0];
    ibiff->value->dest->value->y = pos[1];
    ibiff->value->dest->value->z = pos[2];

    iconn->value->src->value->x = pos[0];
    iconn->value->src->value->y = pos[1];
    iconn->value->src->value->z = pos[2];

    inew->value->src->value->x = pos[0];
    inew->value->src->value->y = pos[1];
    inew->value->src->value->z = pos[2];
}

void save_original_bifurcation_position (struct segment_node *ibiff, double ori_pos[])
{
    ori_pos[0] = ibiff->value->dest->value->x;
    ori_pos[1] = ibiff->value->dest->value->y;
    ori_pos[2] = ibiff->value->dest->value->z;
}

void initialize_best_position_as_middle_point(double best_pos[], const double ori_pos[])
{
    best_pos[0] = ori_pos[0];
    best_pos[1] = ori_pos[1];
    best_pos[2] = ori_pos[2];
}