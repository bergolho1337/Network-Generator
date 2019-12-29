#include "config.h"

User_Data::User_Data ()
{
    this->param = new std::map<std::string,double>();
}

User_Data::~User_Data ()
{
    delete this->function_name;
    delete this->library_name;
    delete this->param;
}

void User_Data::read (const char filename[])
{
    std::ifstream input(filename);
    if (!input)
        fprintf(stderr,"[config] ERROR! Opening user input file '%s'\n",filename);
    else
        fprintf(stdout,"[config] User input file '%s' opened with sucess\n",filename);

    std::string str1, str2, str3;
    while (input >> str1 && str1 != "[main]");

    input >> str1 >> str2 >> this->num_points;

    input >> str1 >> str2 >> str3;
    this->function_name = new std::string(str3);
    input >> str1 >> str2 >> str3;
    this->library_name = new std::string(str3);

    double value;
    while (input >> str1 >> str2 >> value)
    {
        std::string key(str1);

        this->param->insert(std::pair<std::string,double>(key,value));
    }

    input.close();
}

void User_Data::print ()
{
    fprintf(stdout,"%s\n",PRINT_LINE);
    fprintf(stdout,"Number of points = %u\n",this->num_points);
    fprintf(stdout,"Function name = '%s'\n",this->function_name->c_str());
    fprintf(stdout,"Library name = '%s'\n",this->library_name->c_str());
    if (!this->param->empty())
    {
        fprintf(stdout,"%s\n",PRINT_LINE_2);
        fprintf(stdout,"Additional parameters:\n");
        for (auto it = this->param->begin(); it != this->param->end(); ++it)
            printf("\t%s = %g\n",it->first.c_str(),it->second);
    }
    fprintf(stdout,"%s\n",PRINT_LINE);
}