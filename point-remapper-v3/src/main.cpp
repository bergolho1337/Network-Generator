// Author: Lucas Berg
// Program that receives a cloud of points and remap them using a root point as a reference.
// In this version the remapping is done directly in this code without using the "points-to-faces-mapping" program
// First of all, the STL file with the mesh that will have the points remapped is loaded in memory with its points and faces
// Secondly, the link between the points and faces is done by only using the points and faces from the STL without building a map, as it was done in the previous version of the program
// With the links configured, the next step is to build the mapping graph, which will store the indexes of the neighbour points from each point in the grid.
// Lastly, a BFS is performed from the root point index given as input argument. The points will be remapped following the search procedure.

#include <iostream>
#include <cstdio>
#include <cstdlib>

#include "stl-reader/reader.h"
#include "graph/graph.h"

using namespace std;

void write_mapping_to_file (vector<Point> points, vector<Face> faces)
{
    // Write the mapping graph to a file
    // FORMAT:  The first line is the number of points
    //          Followed by the coordinates of the points
    //          Next, each line of the file represents a link between a point and face
    //          <num_points>
    //          <x> <y> <z>
    //              .
    //              .
    //              .
    //          <point_i> <face_j>
    //              .
    //              .
    //   
    FILE *file = fopen("outputs/mapping_graph.txt","w+");

    uint32_t num_points = points.size();
    uint32_t num_faces = faces.size();

    fprintf(file,"%u\n",num_points);
    for (uint32_t i = 0; i < num_points; i++)
    {
        uint32_t index;
        double pos[3];

        index = points[i].id;
        pos[0] = points[i].x;
        pos[1] = points[i].y;
        pos[2] = points[i].z;

        fprintf(file,"%lf %lf %lf\n",pos[0],pos[1],pos[2]);
    } 

    for (uint32_t i = 0; i < num_faces; i++)
    {
        uint32_t point_index;
        uint32_t face_index;

        // Get the index from the current face
        face_index = i;

        // Get the indexes from each vertex 
        point_index = faces[i].vertex_index_1;
        fprintf(file,"%u %u\n",point_index,face_index);
        
        point_index = faces[i].vertex_index_2;
        fprintf(file,"%u %u\n",point_index,face_index);
        
        point_index = faces[i].vertex_index_3;
        fprintf(file,"%u %u\n",point_index,face_index);


    }

    fclose(file);
}

void build_mapping (Graph *g, vector<Point> points, vector<Face> faces, vector<Link> &links)
{

    uint32_t num_points = points.size();
    uint32_t num_faces = faces.size();

    // Insert all the mesh points into the graph
    for (uint32_t i = 0; i < num_points; i++)
    {
        uint32_t index;
        double pos[3];

        index = points[i].id;
        pos[0] = points[i].x;
        pos[1] = points[i].y;
        pos[2] = points[i].z;

        g->insert_node_graph(index,pos);
    }

    for (uint32_t i = 0; i < num_faces; i++)
    {
        uint32_t point_index;
        uint32_t face_index;

        // Get the indexes from each vertex 
        uint32_t face_vertex_ids[3];
        face_vertex_ids[0] = faces[i].vertex_index_1;
        face_vertex_ids[1] = faces[i].vertex_index_2;
        face_vertex_ids[2] = faces[i].vertex_index_3;

        // Get the index from the current face
        face_index = i;

        Link link1(face_vertex_ids[0],face_index);
        Link link2(face_vertex_ids[1],face_index);
        Link link3(face_vertex_ids[2],face_index);

        links.push_back(link1);
        links.push_back(link2);
        links.push_back(link3);
        
    }

    printf("[!] Writing TXT mapping graph file ... File will be saved at 'outputs/mapping_graph.txt'\n");
    write_mapping_to_file(points,faces);
}

void link_points_to_faces (Graph *g, vector<Face> faces, vector<Link> links)
{
    uint32_t num_links = links.size();

    for (uint32_t i = 0; i < num_links; i++)
    {
        printf("\t[!] Working on link number (%u out of %u) ...\n",i,num_links);

        uint32_t point_index = links[i].point_index;
        uint32_t face_index = links[i].face_index;

        //cout << "Link " << point_index << " " << face_index << endl;

        // Get the indexes from the points of the current face
        uint32_t face_points_id[3];
        face_points_id[0] = faces[face_index].vertex_index_1;
        face_points_id[1] = faces[face_index].vertex_index_2;
        face_points_id[2] = faces[face_index].vertex_index_3;

        // Now insert the edges linking the nodes on the graph
        g->insert_edge_graph(point_index,face_points_id[0]);
        g->insert_edge_graph(point_index,face_points_id[1]);
        g->insert_edge_graph(point_index,face_points_id[2]);
    }

    g->print();
}

int main (int argc, char *argv[])
{
    if (argc-1 != 2)
    {
        printf("----------------------------------------------------------------------------\n");
        printf("Usage:> %s <mesh_STL_filename> <root_point_index>\n",argv[0]);
        printf("----------------------------------------------------------------------------\n");
        printf("Example:\n");
        printf("\t%s inputs/generated_grid.stl 0\n",argv[0]);
        printf("\t%s inputs/generated_grid.stl 1300\n",argv[0]);
        printf("\t%s inputs/generated_grid_refined.stl 5100\n",argv[0]);
	    printf("----------------------------------------------------------------------------\n");
        exit(EXIT_FAILURE);   
    }

    // Read the STL mesh file and store its faces in a vector
    printf("[!] Reading STL filename '%s' ...\n",argv[1]);
    vector<Point> points;
    vector<Face> faces;
    read_faces(argv[1],points,faces);
    //print_faces(faces);

    // Build the mapping between points and faces
    printf("[!] Building points to faces links ...\n");
    Graph *the_mapping = new Graph();
    vector<Link> links;
    build_mapping(the_mapping,points,faces,links); 

    // Link the neighbour points using the mapping graph
    printf("[!] Building the mapping graph ...\n");
    link_points_to_faces(the_mapping,faces,links);

    // Run a BFS to remap the numbering from the points of the mesh
    printf("[!] Running a BFS over the graph from the point index '%u' ...\n",atoi(argv[2]));
    uint32_t root_index = atoi(argv[2]);
    vector<Node> mapped_points;
    the_mapping->breadth_first_search(mapped_points,root_index);

    printf("[!] Writing mapped points in a PTS file ... File will be saved at 'outputs/mapped_points.pts'\n");
    write_mapped_points_to_pts(mapped_points);

    printf("[!] Writing mapped points in a TXT file ... File will be saved at 'outputs/mapped_points.pts'\n");
    write_mapped_points_to_txt(mapped_points);

    delete the_mapping;

    return 0;
}
