#include "walker_config.h"

struct walker_config* new_walker_config ()
{
    struct walker_config *result = (struct walker_config*)malloc(sizeof(struct walker_config));

    result->handle = NULL;
    result->library_name = NULL;
    result->params = new std::map<std::string,double>();

    return result;
}

void free_walker_config (struct walker_config *the_walker_config)
{
    delete the_walker_config->params;

    if (the_walker_config->library_name)
        free(the_walker_config->library_name);

    if (the_walker_config->handle)
        dlclose(the_walker_config->handle);
    
    free(the_walker_config);
}

void set_walker_functions (struct walker_config *the_walker_config)
{
    assert(the_walker_config->library_name);

    char *library_path = the_walker_config->library_name;

    // Open library
    the_walker_config->handle = dlopen(library_path,RTLD_LAZY);
    if (!the_walker_config->handle) 
    {
        fprintf(stderr,"%s\n",dlerror());
        exit(EXIT_FAILURE);
    }
    else
    {
        fprintf(stdout,"\n[walker_config] Function library \"%s\" opened with sucess\n",library_path);
    }

    // Search for the 'move' function
    char walker_move_function_name[50] = "move";
    the_walker_config->move_function = (set_walker_move_function_fn*)dlsym(the_walker_config->handle,walker_move_function_name);
    if (dlerror() != NULL)  
    {
        fprintf(stderr, "[walker_config] \"%s\" function not found in the provided function library\n",walker_move_function_name);
        exit(EXIT_FAILURE);
    }
    else
    {
        fprintf(stdout, "[walker_config] Using the function \"%s\" as walker move function\n",walker_move_function_name);
    }

    // Search for the 'respawn' function
    char walker_respawn_function_name[50] = "respawn";
    the_walker_config->respawn_function = (set_walker_respawn_function_fn*)dlsym(the_walker_config->handle,walker_respawn_function_name);
    if (dlerror() != NULL)  
    {
        fprintf(stderr, "[walker_config] \"%s\" function not found in the provided function library\n",walker_respawn_function_name);
        exit(EXIT_FAILURE);
    }
    else
    {
        fprintf(stdout, "[walker_config] Using the function \"%s\" as walker respawn function\n",walker_respawn_function_name);
    }

    // Search for the 'draw' function
    char walker_draw_domain_function_name[50] = "draw";
    the_walker_config->draw_domain_function = (set_walker_draw_domain_function_fn*)dlsym(the_walker_config->handle,walker_draw_domain_function_name);
    if (dlerror() != NULL)  
    {
        fprintf(stderr, "[walker_config] \"%s\" function not found in the provided function library\n",walker_draw_domain_function_name);
        exit(EXIT_FAILURE);
    }
    else
    {
        fprintf(stdout, "[walker_config] Using the function \"%s\" as walker draw domain function\n",walker_draw_domain_function_name);
    }
}

void print_walker_config (struct walker_config *the_walker_config)
{
    printf("Walker function library name = \"%s\"\n",the_walker_config->library_name);

    if (!the_walker_config->params->empty())
    {
        printf("Parameters:\n");
        for (auto it = the_walker_config->params->begin(); it != the_walker_config->params->end(); ++it)
            printf("\t%s = %lf\n",it->first.c_str(),it->second);
    }
}
