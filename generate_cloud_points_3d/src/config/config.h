#ifndef CONFIG_H
#define CONFIG_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <dlfcn.h>

#include <string>
#include <map>
#include <fstream>

#include "../utils/utils.h"

class User_Data
{
public:
    u_int32_t num_points;
    std::string *function_name;
    std::string *library_name;
    
    std::map<std::string,double> *param;
public:
    User_Data ();
    ~User_Data ();
    void read (const char filename[]);
    void print ();
};

#endif