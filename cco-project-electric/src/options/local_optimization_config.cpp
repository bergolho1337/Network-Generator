#include "local_optimization_config.h"

struct local_optimization_config* new_local_optimization_config ()
{
    struct local_optimization_config *result = (struct local_optimization_config*)malloc(sizeof(struct local_optimization_config));

    result->handle = NULL;
    result->name = NULL;
    result->function = NULL;

    result->first_call = true;
    
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
    assert(config);

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

void save_original_bifurcation_position (struct segment_node *iconn, double ori_pos[])
{
    ori_pos[0] = iconn->value->src->value->x;
    ori_pos[1] = iconn->value->src->value->y;
    ori_pos[2] = iconn->value->src->value->z;
}

void initialize_best_position_as_middle_point(double best_pos[], const double ori_pos[])
{
    best_pos[0] = ori_pos[0];
    best_pos[1] = ori_pos[1];
    best_pos[2] = ori_pos[2];
}