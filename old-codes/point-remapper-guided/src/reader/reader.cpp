#include "reader.h"

Reader::Reader (int argc, char *argv[])
{
    char *network_filename = argv[1];
    this->the_network = new Graph(network_filename);
    this->the_network->print();
    //this->the_network->print_terminals();

    char *cloud_points_filename = argv[2];
    this->the_cloud = new Cloud_Point(cloud_points_filename);
    //this->the_cloud->print();

}

void Reader::remap_points_using_graph ()
{
    uint32_t total_nodes = this->the_network->get_total_nodes();
    int *parents = new int[total_nodes]();

    int source_index = 0;
    //int source_index = 142;
    Node *source_node = this->the_network->search_node(source_index);
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

        check_cloud_points_inside_node(u);

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
// DFS
/*
    depth_first_search(source_index);

    // Only the network pathways
    write_remapped_points_to_vtk();
    write_points_to_pts();

    for (uint32_t i = 0; i < this->the_cloud->the_points.size(); i++)
    {
        bool is_taken = this->the_cloud->taken[i];

        if (!is_taken)
        {
            Point_3d point = this->the_cloud->the_points[i];   

            this->the_cloud->the_remapped_points.push_back(point);
        }
    }
*/

    // All cloud of points remapped
    write_remapped_points_to_vtk();
    write_points_to_pts();

}

void Reader::check_cloud_points_inside_node (Node *u)
{
    double center[3];
    center[0] = u->x;
    center[1] = u->y;
    center[2] = u->z;

    uint32_t num_points = this->the_cloud->the_points.size();
    for (uint32_t i = 0; i < num_points; i++)
    {
        Point_3d point = this->the_cloud->the_points[i];
        bool is_taken = this->the_cloud->taken[i];

        if (is_inside(point,center) && !is_taken)
        {
            this->the_cloud->the_remapped_points.push_back(point);

            this->the_cloud->taken[i] = true;
        }
    }
}

bool Reader::is_inside (Point_3d point, const double center[])
{
    double r = RADIUS;
    double dx = (point.x - center[0]);
    double dy = (point.y - center[1]);
    double dz = (point.z - center[2]);

    double value = ((dx*dx)) + ((dy*dy)) + ((dz*dz));
    if (value < r*r)
        return true;
    else
        return false;
}

void Reader::write_remapped_points_to_vtk ()
{
    char filename[200];
    sprintf(filename,"outputs/remapped_points.vtk");
    FILE *file = fopen(filename,"w+");

    fprintf(file,"# vtk DataFile Version 3.0\n");
    fprintf(file,"Cloud\n");
    fprintf(file,"ASCII\n");
    fprintf(file,"DATASET POLYDATA\n");
    fprintf(file,"POINTS %lu float\n",this->the_cloud->the_remapped_points.size());
    for (uint32_t i = 0; i < this->the_cloud->the_remapped_points.size(); i++)
        fprintf(file,"%lf %lf %lf\n",this->the_cloud->the_remapped_points[i].x*SCALE_RATIO,this->the_cloud->the_remapped_points[i].y*SCALE_RATIO,this->the_cloud->the_remapped_points[i].z*SCALE_RATIO);
        //fprintf(file,"%lf %lf %lf\n",this->the_cloud->the_remapped_points[i].x,this->the_cloud->the_remapped_points[i].y,this->the_cloud->the_remapped_points[i].z);
    fprintf(file,"VERTICES %lu %lu\n",this->the_cloud->the_remapped_points.size(),this->the_cloud->the_remapped_points.size()*2);
    for (uint32_t i = 0; i < this->the_cloud->the_remapped_points.size(); i++)
        fprintf(file,"1 %u\n",i);
    
    fclose(file);
}

void Reader::write_points_to_pts ()
{
    char filename[200];
    sprintf(filename,"outputs/remapped_points.pts");
    FILE *file = fopen(filename,"w+");

    fprintf(file,"%lu\n",this->the_cloud->the_remapped_points.size());
    for (uint32_t i = 0; i < this->the_cloud->the_remapped_points.size(); i++)
        fprintf(file,"%lf %lf %lf\n",this->the_cloud->the_remapped_points[i].x*SCALE_RATIO,this->the_cloud->the_remapped_points[i].y*SCALE_RATIO,this->the_cloud->the_remapped_points[i].z*SCALE_RATIO);
        //fprintf(file,"%lf %lf %lf\n",this->the_cloud->the_remapped_points[i].x*0.0005,this->the_cloud->the_remapped_points[i].y*0.0005,this->the_cloud->the_remapped_points[i].z*0.0005);
    
    fclose(file);
}


void Reader::depth_first_search (const uint32_t source_index)
{
    Node *source_node = this->the_network->search_node(source_index);
    uint32_t total_nodes = this->the_network->get_total_nodes();

    vector<int> dfs_num;
    dfs_num.assign(total_nodes,DFS_WHITE);

    dfs(source_node,dfs_num);
}

void Reader::dfs (Node *u, vector<int> &dfs_num)
{
    int u_index = u->index;
    Edge *v = u->list_edges;

    check_cloud_points_inside_node(u);

    dfs_num[u_index] = DFS_BLACK;
    //printf("[depth_first_search] Current on index %d\n",u_index);

    while (v != NULL)
    {
        int v_index = v->index;
        if (dfs_num[v_index] == DFS_WHITE)
        {
            dfs(v->dest,dfs_num);
        }
        v = v->next;
    }
}