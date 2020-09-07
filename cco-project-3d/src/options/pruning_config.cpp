#include "pruning_config.h"

struct pruning_config* new_pruning_config ()
{
    struct pruning_config *result = (struct pruning_config*)malloc(sizeof(struct pruning_config));

    result->handle = NULL;
    result->function_name = NULL;
    result->library_name = NULL;
    result->function = NULL;
    
    result->params = new std::map<std::string,double>();

    return result;
}

void free_pruning_config (struct pruning_config *config)
{
    //delete config->params;
    if (config->library_name)
        free(config->library_name);

    if (config->function_name)
        free(config->function_name);

    if (config->handle)
        dlclose(config->handle);
    
    free(config);
}

void set_pruning_function (struct pruning_config *config)
{
    assert(config);

    char library_path[200];
    if (config->library_name)
        strcpy(library_path,config->library_name);
    else
        strcpy(library_path,"./shared-libs/libdefault_pruning.so");

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
    
    char *pruning_function_name = config->function_name;

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
    printf("Pruning function name = \"%s\"\n",config->function_name);

    if (!config->params->empty())
    {
        printf("Parameters:\n");
        for (auto it = config->params->begin(); it != config->params->end(); ++it)
            printf("\t%s = %lf\n",it->first.c_str(),it->second);
    }
}
