#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>

#include <map>
#include <vector>

#include "stl-reader/reader.h"
#include "graph/graph.h"
#include "utils/utils.h"

const double RADIUS = 0.001;

bool is_inside (Point p, const double center[])
{
    double r = RADIUS;
    double dx = (p.x - center[0]);
    double dy = (p.y - center[1]);
    double dz = (p.z - center[2]);

    double value = ((dx*dx)) + ((dy*dy)) + ((dz*dz));
    if (value < r*r)
        return true;
    else
        return false;
}

void check_faces_inside_node (Node *u, vector<Point> points, vector<Face> &faces)
{
    double center[3];
    center[0] = u->x;
    center[1] = u->y;
    center[2] = u->z;

    uint32_t num_faces = faces.size();
    for (uint32_t i = 0; i < num_faces; i++)
    {
        Point p1 = points[faces[i].vertex_index_1];
        Point p2 = points[faces[i].vertex_index_2];
        Point p3 = points[faces[i].vertex_index_3];

        if (is_inside(p1,center) || is_inside(p2,center) || is_inside(p3,center))
        {
            faces[i].is_taken = true;
        }
    }
}

void remap_mesh_using_network (vector<Point> points, vector<Face> &faces, Graph *network, const char filename[])
{
    uint32_t total_nodes = network->get_total_nodes();
    int *parents = new int[total_nodes]();

    int source_index = 0;
    Node *source_node = network->search_node(source_index);
    parents[source_index] = -1;

    map<int,int> dist;              // Distance from source to the other nodes
    dist[source_index] = 0;         // Distance from source to source is zero

// BFS
    queue<Node*> q;
    q.push(source_node);            // Enqueue the source first

    while (!q.empty())
    {
        Node *u = q.front(); q.pop();
        int u_index = u->index;

        check_faces_inside_node(u,points,faces);

        Edge *v = u->list_edges;
        while (v != NULL)
        {
            int v_index = v->index;
            if (!dist.count(v_index))
            {
                dist[v_index] = dist[u_index] + 1;
                parents[v_index] = u_index;
                q.push(v->dest);
            }
            v = v->next;
        }
    }

    write_faces_to_stl(points,faces,filename);
}

int main (int argc, char *argv[])
{
    if (argc-1 != 3)
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
    char *network_filename = argv[2];
    if (!check_file_extension(network_filename,"vtk"))
    {
        fprintf(stderr,"[-] ERROR! Invalid network format! The network must be in the VTK format.\n");
        exit(EXIT_FAILURE);
    }
    char *output_filename = argv[3];
    if (!check_file_extension(output_filename,"stl"))
    {
        fprintf(stderr,"[-] ERROR! Invalid mesh format! The mesh must be in the STL format.\n");
        exit(EXIT_FAILURE);
    }

    vector<Point> points;
    vector<Face> faces;
    read_faces_from_stl(mesh_filename,points,faces);

    Graph *network = new Graph(network_filename);

    remap_mesh_using_network(points,faces,network,output_filename);
    
    delete network;

    return 0;
}
