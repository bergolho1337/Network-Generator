#include "utils.h"
                              
bool check_file_extension (const char filename[], const char extension[])
{
  uint32_t i, j;
  uint32_t size = strlen(filename);

  char file_extension[5];
  for (i = size-3, j = 0; i < size; i++)
  {
    file_extension[j] = filename[i];
    j++;
  }
  file_extension[j] = '\0';

  if (strcmp(file_extension,extension) == 0)
    return true;
  else
    return false;
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
  printf("Usage:> %s <input_filename> <output_filename>\n",pname);
  printf("%s\n",PRINT_LINE);
}