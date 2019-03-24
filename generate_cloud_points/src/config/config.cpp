#include "config.h"

struct user_data* new_user_data ()
{
    struct user_data *config = (struct user_data*)malloc(sizeof(struct user_data));

    config->param = new std::map<std::string,double>();

    return config; 
}

void free_user_data (struct user_data *config)
{
    delete config->function_name;
    delete config->library_name;

    delete config->param;

    free(config);
}

void read_user_input (struct user_data *config, const char filename[])
{
    
    std::ifstream input(filename);
    if (!input)
        fprintf(stderr,"[config] ERROR! Opening user input file '%s'\n",filename);
    else
        fprintf(stdout,"[config] User input file '%s' opened with sucess\n",filename);

    std::string str1, str2, str3;
    while (input >> str1 && str1 != "[main]");

    input >> str1 >> str2 >> config->area;
    input >> str1 >> str2 >> config->num_points;

    input >> str1 >> str2 >> str3;
    config->function_name = new std::string(str3);
    input >> str1 >> str2 >> str3;
    config->library_name = new std::string(str3);

    double value;
    while (input >> str1 >> str2 >> value)
    {
        std::string key(str1);

        config->param->insert(std::pair<std::string,double>(key,value));
    }

    input.close();
}

void print_user_input (struct user_data *config)
{
    fprintf(stdout,"%s\n",PRINT_LINE);
    fprintf(stdout,"Area = %g\n",config->area);
    fprintf(stdout,"Number of points = %u\n",config->num_points);
    fprintf(stdout,"Function name = '%s'\n",config->function_name->c_str());
    fprintf(stdout,"Library name = '%s'\n",config->library_name->c_str());
    if (!config->param->empty())
    {
        fprintf(stdout,"%s\n",PRINT_LINE_2);
        fprintf(stdout,"Additional parameters:\n");
        for (auto it = config->param->begin(); it != config->param->end(); ++it)
            printf("\t%s = %g\n",it->first.c_str(),it->second);
    }
    fprintf(stdout,"%s\n",PRINT_LINE);
}