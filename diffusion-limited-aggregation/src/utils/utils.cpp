#include "utils.h"

double calculate_distance (const double p1[], const double p2[])
{
    return (p2[0] - p1[0])*(p2[0] - p1[0]) + (p2[1] - p1[1])*(p2[1] - p1[1]) + (p2[2] - p1[2])*(p2[2] - p1[2]);
}

double calculate_euclidean_norm (const double x1, const double y1, const double z1,\
                                const double x2, const double y2, const double z2)
{
    return sqrt( pow(x2-x1,2) + pow(y2-y1,2) + pow(z2-z1,2) ); 
}

double generate_random_number ()
{
    uint8_t sign = rand() % 2;
    double value = (double) rand() / (double) RAND_MAX;
    if (sign)
        return value;
    else
        return -value;
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
        //fprintf(stderr,"Not found \"%s\" ... using default value\n",key.c_str());
        return false;
    }
}
                              

void print_progress_bar (const uint32_t cur_iter, const uint32_t max_iter)
{
  double percentage = (cur_iter / (double) max_iter) * 100.0;
  uint32_t filled_length = nearbyint(100 * cur_iter / max_iter);
  
  std::string the_bar;
  for (uint32_t i = 0; i < filled_length; i++)
    the_bar += "\u2588";    // Using Unicode here ... (filled box)
  for (uint32_t i = 0; i < 100-filled_length; i++)
    the_bar += "-";
  
  printf("%s |%s| %.1lf%% %s\r","Progress",the_bar.c_str(),percentage,"Complete"); // Carriage return
  fflush(stdout); // Flush the standard output
}

void usage (const char pname[])
{
  printf("%s\n",PRINT_LINE);
  printf("Usage:> %s <input_file>\n",pname);
  printf("%s\n",PRINT_LINE);
}