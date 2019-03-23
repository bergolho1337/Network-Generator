#include "utils.h"

void usage (const char pname[])
{
    fprintf(stderr,"%s\n",PRINT_LINE);
    fprintf(stderr,"Usage:> %s <input_config_file>");
    fprintf(stderr,"%s\n",PRINT_LINE_2);
    fprintf(stderr,"<input_config_file> = Input file with user configurations\n");
    fprintf(stderr,"%s\n",PRINT_LINE);
}