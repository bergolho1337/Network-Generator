#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>

#include <map>
#include <vector>

#include "stl-reader/reader.h"
#include "utils/utils.h"

void build_mapping (vector< vector<uint32_t> > &the_mapping, vector<Point> points, vector<Face> faces)
{
    uint32_t num_points = points.size();
    uint32_t num_faces = faces.size();

    // Alocate the nodes
    the_mapping.assign(num_points,vector<uint32_t>());

    // Pass through each face
    for (uint32_t i = 0; i < num_faces; i++)
    {
        uint32_t face_index = i;

        uint32_t vertex_index_1 = faces[i].vertex_index_1;
        uint32_t vertex_index_2 = faces[i].vertex_index_2;
        uint32_t vertex_index_3 = faces[i].vertex_index_3;

        // Add the edge: point_index -> face_index
        the_mapping[vertex_index_1].push_back(face_index);
        the_mapping[vertex_index_2].push_back(face_index);
        the_mapping[vertex_index_3].push_back(face_index);
    }
}

// Write the mapping graph to a file
    // FORMAT:  The first line is the number of points
    //          Next, each line of the file represents a link between a point and face
    //          <num_points>
    //          <point_i> <face_j>
    //              .
    //              .
    //              .
void write_mapping_to_txt (vector< vector<uint32_t> > the_mapping)
{
    FILE *file = fopen("outputs/mapping_graph.txt","w+");

    fprintf(file,"%u\n",the_mapping.size());
    for (uint32_t i = 0; i < the_mapping.size(); i++)
    {
        for (uint32_t j = 0; j < the_mapping[i].size(); j++)
        {
            uint32_t point_index = i;
            uint32_t face_index = the_mapping[i][j];

            fprintf(file,"%u %u\n",point_index,face_index);
        }
    }

    fclose(file);
}

int main (int argc, char *argv[])
{
    if (argc-1 != 1)
    {
        usage(argv[0]);
        exit(EXIT_FAILURE);
    }

    // Get the input mesh file in STL
    char *mesh_filename = argv[1];
    if (!check_file_extension(mesh_filename,"stl"))
    {
        fprintf(stderr,"[-] ERROR! Invalid mesh format! The mesh must be in the STL format.\n");
        exit(EXIT_FAILURE);
    }

    vector<Point> points;
    vector<Face> faces;
    read_faces_from_stl(mesh_filename,points,faces);

    vector< vector<uint32_t> > the_mapping;
    build_mapping(the_mapping,points,faces);

    write_mapping_to_txt(the_mapping);

    return 0;
}
