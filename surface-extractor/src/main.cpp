// Author: Lucas Berg
// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// Program that given a mesh in STL format reads its points and faces, followed by a mesh graph in TXT and a PTS containing the remapped points given by the Point-Remapper-v3,
// returns the surface which have the points from [0,id_limit] of the PTS file.
// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#include <iostream>
#include <cstdio>
#include <cstdlib>

#include "stl-reader/reader.h"

using namespace std;

void mark_faces_within_limit (vector< vector<uint32_t> > the_mapping, vector<Point> mapped_points, const uint32_t index_limit, vector<Face> &faces)
{
    // Pass through the remapped points vector
    for (uint32_t i = 0; i < mapped_points.size(); i++)
    {
        // Only check the points within the interval [0,id_limit]
        if (i <= index_limit)
        {
            // Get the original index from the point in the STL file ...
            uint32_t index = mapped_points[i].id;

            // Now using the mesh graph, check the list of faces and mark them to be written in the output STL.
            for (uint32_t j = 0; j < the_mapping[index].size(); j++)
            {
                uint32_t face_index = the_mapping[index][j];
  
                faces[face_index].marked = true;
            }
        }
    }
}

void read_mapping_from_txt(vector< vector<uint32_t> > &the_mapping, const char filename[])
{
    FILE *file = fopen(filename,"r");

    uint32_t num_nodes;
    fscanf(file,"%u",&num_nodes);
    the_mapping.assign(num_nodes,vector<uint32_t>());

    for (uint32_t i = 0; i < num_nodes; i++)
    {
        double pos[3];
        fscanf(file,"%lf %lf %lf",&pos[0],&pos[1],&pos[2]);
    }

    uint32_t link[2];
    // Link: id_point -> id_face
    while (fscanf(file,"%u %u",&link[0],&link[1]) != EOF)
    {
        //printf("%u %u\n",link[0],link[1]);
        the_mapping[link[0]].push_back(link[1]);
    }

    fclose(file);
}

// Remember that the numeration from the points in the "mapped_points.pts" is different from the original STL ...
void read_mapped_points_from_txt (vector<Point> &mapped_points, const char filename[])
{
    FILE *file = fopen(filename,"r");

    uint32_t num_points;
    fscanf(file,"%u",&num_points);

    uint32_t index;
    double pos[3];
    for (uint32_t i = 0; i < num_points; i++)
    {
        // Reading format for the points:
        //      "original_index_from_stl" "x" "y" "z"
        fscanf(file,"%u %lf %lf %lf",&index,&pos[0],&pos[1],&pos[2]);

        Point p(index,pos[0],pos[1],pos[2]);

        mapped_points.push_back(p);
    }

    fclose(file);
}

int main (int argc, char *argv[])
{
    if (argc-1 != 4)
    {
        printf("------------------------------------------------------------------------------------------------------------------------------------------\n");
        printf("Usage:> %s <mesh_STL_filename> <mesh_graph_TXT_filename> <mapped_points_with_ids_TXT_filename> <index_limit>\n",argv[0]);
        printf("------------------------------------------------------------------------------------------------------------------------------------------\n");
        printf("Example:\n");
        printf("\t%s inputs/elizabeth_biventricular_full.stl inputs/elizabeth_mapping_graph_for_LV.txt inputs/elizabeth_mapped_points_with_id_for_LV.txt 150000\n",argv[0]);
        printf("\t%s inputs/elizabeth_biventricular_full.stl inputs/elizabeth_mapping_graph_for_RV.txt inputs/elizabeth_mapped_points_with_id_for_RV.txt 200000\n",argv[0]);
        printf("------------------------------------------------------------------------------------------------------------------------------------------\n");
        exit(EXIT_FAILURE);   
    }

    // Read the STL mesh file and store its faces in a vector
    printf("[!] Reading STL filename '%s' ...\n",argv[1]);
    vector<Point> points;
    vector<Face> faces;
    read_faces(argv[1],points,faces);
    //print_faces(faces);

    // Read the TXT mapping graph
    printf("[!] Reading TXT with the mapping graph ...\n");
    vector< vector<uint32_t> > the_mapping;
    read_mapping_from_txt(the_mapping,argv[2]);

    // Read
    printf("[!] Reading TXT with the mapped points ...\n");
    vector<Point> mapped_points;
    read_mapped_points_from_txt(mapped_points,argv[3]);

    // Read the index of the maximum point that we should read
    printf("[!] Reading index limit ...\n");
    uint32_t index_limit = atoi(argv[4]);

    // Mark the faces which have the points inside the interval: [0,index_limit]
    printf("[!] Marking faces ...\n");
    mark_faces_within_limit(the_mapping,mapped_points,index_limit,faces);

    printf("[!] Writing output STL file ... File will be saved at 'outputs/processed_mesh.stl'\n");
    write_marked_faces_to_stl(points,faces);

    return 0;
}
