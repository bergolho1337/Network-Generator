#include "utils.h"

void usage (const char pname[])
{
    printf("%s\n",PRINT_LINE);
    printf("Usage:> %s <input_file>\n",pname);
    printf("\t<input_file> = Input file with the graph configuration\n");
    printf("%s\n",PRINT_LINE);
}

double calc_norm (const double x1, const double y1, const double z1,\
                  const double x2, const double y2, const double z2)
{
    return sqrt(pow(x2-x1,2) + pow(y2-y1,2) + pow(z2-z1,2));
}

void get_file_extension (const char filename[], char extension_name[])
{
    uint32_t size = strlen(filename);
    for (uint32_t i = size-3, j = 0; i < size; i++, j++)
        extension_name[j] = filename[i];
    extension_name[size-1] = '\0';
}

bool check_file_extension (const char filename[], const char extension_name[])
{
    char file_extension[4];
    get_file_extension(filename,file_extension);
    
    if (strcmp(file_extension,extension_name) == 0)
        return true;
    else
        return false;
}

double calc_angle_between_vectors (const double u[], const double v[])
{
    double dot_product = u[0]*v[0] + u[1]*v[1] + u[2]*v[2];

    double angle_radians = acos(dot_product);

    // Return the angle in degrees
    return angle_radians * 180.0 / M_PI;
}