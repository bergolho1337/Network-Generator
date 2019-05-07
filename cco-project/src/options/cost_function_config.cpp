#include "cost_function_config.h"

struct cost_function_config* new_cost_function_config ()
{
    struct cost_function_config *result = (struct cost_function_config*)malloc(sizeof(struct cost_function_config));

    result->params = new std::map<std::string,double>();

    return result;
}

void free_cost_function_config (struct cost_function_config *config)
{
    delete config->params;
    free(config->function_name);
    free(config->library_name);

    if (config->handle)
        dlclose(config->handle);
    
    free(config);
}

bool get_parameter_value_from_map (std::map<std::string,double> *params,\
                                const std::string key, double *value)
{
    auto it = params->find(key);

    if (it != params->end())
    {
        //printf("Found \"%s\" on parameter list: %g\n",key.c_str(),it->second);
        *value = it->second;
        return true;
    }
    else
    {
        fprintf(stderr,"Not found \"%s\" \n",key.c_str());
        return false;
    }
}

void set_cost_function (struct cost_function_config *config)
{
    char *library_path = config->library_name;

    config->handle = dlopen(library_path,RTLD_LAZY);
    if (!config->handle) 
    {
        fprintf(stderr,"%s\n",dlerror());
        exit(EXIT_FAILURE);
    }
    else
    {
        fprintf(stdout,"\n[cost_function] Cost function library \"%s\" opened with sucess\n",library_path);
    }
    
    char *cost_function_name = config->function_name;

    config->function = (set_cost_function_fn*)dlsym(config->handle,cost_function_name);
    if (dlerror() != NULL)  
    {
        fprintf(stderr, "[cost_function] \"%s\" function not found in the provided cost function library\n",cost_function_name);
        exit(EXIT_FAILURE);
    }
    else
    {
        fprintf(stdout, "[cost_function] Using the function \"%s\" as cost function\n",cost_function_name);
    }
}

void print_cost_function_config (struct cost_function_config *config)
{
    printf("Cost function library name = \"%s\"\n",config->library_name);
    printf("Cost function name = \"%s\"\n",config->function_name);

    if (!config->params->empty())
    {
        printf("Parameters:\n");
        for (auto it = config->params->begin(); it != config->params->end(); ++it)
            printf("\t%s = %lf\n",it->first.c_str(),it->second);
    }
}