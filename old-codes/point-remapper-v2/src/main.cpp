// Author: Lucas Berg
// This second version of the Point-Remapper program performs the remapping using a graph instead of a growing volume.
// The principle is to first read the mesh 

// Program that receives a cloud of points and remap them using a root point as a reference.
// The remapping is done by envoloping the root point with a sphere, which at each iteration will
// increase its radius. By using this idea, all the current points of the cloud that are inside the spehre and 
// that are not visited yet will be placed in a array of points that can later be writen in .pts format.

#include <iostream>
#include <cstdio>
#include <cstdlib>

#include "stl-reader/reader.h"
#include "graph/graph.h"

using namespace std;

void read_graph_from_stl_and_mapping (Graph *g, vector<Face> faces, map<Point_Table,uint32_t> point_table, const char mapping_filename[])
{
    FILE *file = fopen(mapping_filename,"r");
    if (!file)
    {
        fprintf(stderr,"[-] ERROR! Cannot open filename '%s'\n",mapping_filename);
        exit(EXIT_FAILURE);
    }

    uint32_t num_nodes;
    fscanf(file,"%u",&num_nodes);

    // Insert all the unique points from the mesh
    double pos[3];   
    uint32_t point_id, face_id; 
    map<Point_Table,uint32_t>::iterator it;
    for (it = point_table.begin(); it != point_table.end(); ++it)
    {
        pos[0] = it->first.x;
        pos[1] = it->first.y;
        pos[2] = it->first.z;
        point_id = it->second;
        
        g->insert_node_graph(point_id,pos);
    }

    // Read all the point to face links
    while (fscanf(file,"%u %u",&point_id,&face_id) != EOF)
    {
        // Get the indexes from the points of the current face
        uint32_t face_points_id[3];
        face_points_id[0] = faces[face_id].v1->id;
        face_points_id[1] = faces[face_id].v2->id;
        face_points_id[2] = faces[face_id].v3->id;

        // Now insert the edges linking the nodes on the graph
        g->insert_edge_graph(point_id,face_points_id[0]);
        g->insert_edge_graph(point_id,face_points_id[1]);
        g->insert_edge_graph(point_id,face_points_id[2]);
    }

    
    //g->print();
    fclose(file);
}

int main (int argc, char *argv[])
{
    if (argc-1 != 3)
    {
        printf("----------------------------------------------------------------------------\n");
        printf("Usage:> %s <mesh_STL_filename> <mapping_TXT_filename> <root_point_index>\n",argv[0]);
        printf("----------------------------------------------------------------------------\n");
        printf("Example:\n");
        printf("\t%s inputs/generated_grid.stl inputs/generated_grid_mapping.txt 0\n",argv[0]);
        printf("\t%s inputs/generated_grid.stl inputs/generated_grid_mapping.txt 1300\n",argv[0]);
        printf("----------------------------------------------------------------------------\n");
        exit(EXIT_FAILURE);   
    }

    // Read the STL mesh file and store its faces in a vector
    vector<Face> faces;
    map<Point_Table,uint32_t> point_table;
    printf("[!] Reading STL mesh file and building points table ...\n");
    read_faces(argv[1],faces,point_table);
    //print_faces(faces);

    // Read the mapping file and build the graph
    Graph *the_mapping = new Graph();
    printf("[!] Reading TXT graph mapping file ...\n");
    read_graph_from_stl_and_mapping(the_mapping,faces,point_table,argv[2]);

    // Run a BFS to remap the numbering from the points of the mesh
    uint32_t root_index = atoi(argv[3]);
    vector<Node> mapped_points;
    printf("[!] Running BFS over the graph ...\n");
    the_mapping->breadth_first_search(mapped_points,root_index);

    printf("[!] Writing PTS mapping output file ...\n");
    write_mapped_points_to_pts(mapped_points);

    delete the_mapping;

    return 0;
}
