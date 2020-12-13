#include "reader/reader.h"
#include "graph/graph.h"

int main (int argc, char *argv[])
{
    if (argc-1 != 1)
    {
        printf("===============================================================================================\n");
        printf("Usage:> %s <input_purkinje_network>\n",argv[0]);
        printf("===============================================================================================\n");
        exit(EXIT_FAILURE);
    }

    std::string purkinje_filename = argv[1];

    // [PURKINJE]
    VTK_Reader *reader = new VTK_Reader(purkinje_filename);

    Graph *network = new Graph(reader->the_points,reader->the_lines);    

    network->write_info();

    return 0;
}