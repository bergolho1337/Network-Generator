// This program takes two cloud of points and merge them

#include "reader/reader.h"

int main (int argc, char *argv[])
{
    if (argc-1 != 3)
    {
        printf("===============================================================================================\n");
        printf("Usage:> %s <input_cloud_1> <input_cloud_2> <output_cloud>\n",argv[0]);
        printf("===============================================================================================\n");
        printf("<input_cloud_1> = Input cloud representing the inactives points\n");
        printf("<input_cloud_2> = Input cloud representing the actives points\n");
        printf("<output_cloud> = Output merged cloud\n");
        printf("===============================================================================================\n");
        exit(EXIT_FAILURE);
    }

    std::string input_cloud_filename_1 = argv[1];
    std::string input_cloud_filename_2 = argv[2];
    std::string output_filename = argv[3];

    VTK_Reader *inactive_cloud = new VTK_Reader(input_cloud_filename_1);
    VTK_Reader *active_cloud = new VTK_Reader(input_cloud_filename_2);

    inactive_cloud->concatenate(active_cloud);

    inactive_cloud->write(output_filename);

    return 0;
}
