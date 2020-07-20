#include "pruning_config.h"

struct pruning_config* new_pruning_config ()
{
    struct pruning_config *result = (struct pruning_config*)malloc(sizeof(struct pruning_config));

    result->handle = NULL;
    result->name = NULL;
    result->function = NULL;
    
    result->params = new std::map<std::string,double>();

    return result;
}

void free_pruning_config (struct pruning_config *config)
{
    //delete config->params;
    if (config->name)
        free(config->name);

    if (config->handle)
        dlclose(config->handle);
    
    free(config);
}

void set_pruning_function (struct pruning_config *config)
{
    assert(config);

    char library_path[200] = "./shared-libs/libdefault_pruning.so";

    config->handle = dlopen(library_path,RTLD_LAZY);
    if (!config->handle) 
    {
        fprintf(stderr,"%s\n",dlerror());
        exit(EXIT_FAILURE);
    }
    else
    {
        fprintf(stdout,"\n[pruning] Pruning library \"%s\" opened with sucess\n",library_path);
    }
    
    char *pruning_function_name = config->name;

    config->function = (set_pruning_function_fn*)dlsym(config->handle,pruning_function_name);
    if (dlerror() != NULL)  
    {
        fprintf(stderr, "[pruning] \"%s\" function not found in the provided pruning library\n",pruning_function_name);
        exit(EXIT_FAILURE);
    }
    else
    {
        fprintf(stdout, "[pruning] Using the function \"%s\" as pruning function\n",pruning_function_name);
    }
}

void print_pruning_function_config (struct pruning_config *config)
{
    printf("Pruning function name = \"%s\"\n",config->name);

    if (!config->params->empty())
    {
        printf("Parameters:\n");
        for (auto it = config->params->begin(); it != config->params->end(); ++it)
            printf("\t%s = %lf\n",it->first.c_str(),it->second);
    }
}
