
#ifndef UTILS_H
#define UTILS_H

#include <cstdio>
#include <cstdlib>
#include <cstdbool>
#include <cstdint>
#include <cmath>
#include <cstring>

#include <vector>
#include <map>
#include <string>

#define PRINT_LINE "==========================================================================================="
#define PRINT_DOTS "..........................................................................................."

void usage (const char pname[]);
void print_progress_bar (const uint32_t cur_iter, const uint32_t max_iter);
bool check_file_extension (const char filename[], const char extension[]);

#endif