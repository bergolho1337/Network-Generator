// Author: Lucas Berg

#include <iostream>
#include <vector>
#include <string>
#include <algorithm>

#include <cstdio>
#include <cstdlib>
#include <cmath>

using namespace std;

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

int main (int argc, char *argv[])
{
    uint32_t max_iter = 1000000;

    for (uint32_t iter = 0; iter < max_iter; iter++)
    {
      print_progress_bar(iter,max_iter);
    }
    printf("\n");

    return 0;
}
