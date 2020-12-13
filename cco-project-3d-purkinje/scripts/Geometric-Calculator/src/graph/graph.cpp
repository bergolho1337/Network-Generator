#include "graph.h"

Graph::Graph ()
{

}

Graph::Graph (std::vector<Point> points, std::vector<Line> lines)
{
    uint32_t num_nodes = points.size();
    uint32_t num_edges = lines.size();

    this->list_nodes.assign(num_nodes,Node());
    for (uint32_t i = 0; i < num_nodes; i++)
    {
        uint32_t id = points[i].id;
        double x = points[i].x;
        double y = points[i].y;
        double z = points[i].z;

        this->list_nodes[i].id = id;
        this->list_nodes[i].x = x;
        this->list_nodes[i].y = y;
        this->list_nodes[i].z = z;
    }


    for (uint32_t i = 0; i < num_edges; i++)
    {
        uint32_t src = lines[i].src;
        uint32_t dest = lines[i].dest;
        double length = calc_norm(points[src].x,points[src].y,points[src].z,\
                                points[dest].x,points[dest].y,points[dest].z);
        
        Edge e(dest,length);
        this->list_nodes[src].list_edges.push_back(e);
    }

    total_nodes = num_nodes;
    total_edges = num_edges;
}

void Graph::depth_first_search (const uint32_t src_id, std::vector<double> &the_segments)
{
    std::vector<bool> dfs_visited;
    dfs_visited.assign(this->total_nodes,false);

    uint32_t total_segments = 0;
    double segment_size = 0.0;
    uint32_t flag = 0;

    dfs(this->list_nodes[src_id],dfs_visited,the_segments,segment_size,flag,total_segments);
}

void Graph::dfs (Node u, std::vector<bool> &dfs_visited, std::vector<double> &segments, double &segment_size, uint32_t &flag, uint32_t &total_segments)
{
    flag = 0;

    uint32_t u_index = u.id;
    dfs_visited[u_index] = true;

    uint32_t num_edges = u.list_edges.size();
    for (uint32_t j = 0; j < num_edges; j++)
    {
        uint32_t v_index = u.list_edges[j].dest_id;
        double length = u.list_edges[j].length;

        if (!dfs_visited[v_index])
        {
            flag++;
            segment_size += length;
            dfs(this->list_nodes[v_index],dfs_visited,segments,segment_size,flag,total_segments);
        }
    }
    // Terminal edge
    if (u.list_edges.size() == 0)
    {
        total_segments++;
        segments.push_back(segment_size);
        segment_size = 0.0;
    }
    if (flag == 2)
    {
        total_segments++;
    }
}

void Graph::print ()
{
    for (uint32_t i = 0; i < this->total_nodes; i++)
    {
        this->list_nodes[i].print();
    }
}

void Graph::get_edges (std::vector<double> &the_edges)
{
    for (uint32_t i = 0; i < this->list_nodes.size(); i++)
        for (uint32_t j = 0; j < this->list_nodes[i].list_edges.size(); j++)
            the_edges.push_back(this->list_nodes[i].list_edges[j].length*UM_TO_MM);
}

void Graph::get_segments (std::vector<double> &the_segments)
{
    depth_first_search(0,the_segments);
}

void Graph::get_bifurcation_angles (std::vector<double> &the_angles)
{
    for (uint32_t i = 0; i < this->list_nodes.size(); i++)
    {
        if (this->list_nodes[i].list_edges.size() == 2)
        {
            uint32_t src = i;
            uint32_t dest_1 = this->list_nodes[i].list_edges[0].dest_id;
            uint32_t dest_2 = this->list_nodes[i].list_edges[1].dest_id;

            double src_pos[3], dest_pos_1[3], dest_pos_2[3];
            src_pos[0] = this->list_nodes[src].x;
            src_pos[1] = this->list_nodes[src].y;
            src_pos[2] = this->list_nodes[src].z;
            dest_pos_1[0] = this->list_nodes[dest_1].x;
            dest_pos_1[1] = this->list_nodes[dest_1].y;
            dest_pos_1[2] = this->list_nodes[dest_1].z;
            dest_pos_2[0] = this->list_nodes[dest_2].x;
            dest_pos_2[1] = this->list_nodes[dest_2].y;
            dest_pos_2[2] = this->list_nodes[dest_2].z;

            double u1[3], u2[3];
            build_unitary_vector(u1,src_pos,dest_pos_1);
            build_unitary_vector(u2,src_pos,dest_pos_2);

            double angle = calc_angle_between_vectors(u1,u2);
            the_angles.push_back(angle);
        }
    }
}

uint32_t Graph::get_closest_point (Point p)
{
    uint32_t closest_index = 0;
    double closest_dist = __DBL_MAX__;
    for (uint32_t i = 0; i < this->list_nodes.size(); i++)
    {
        double dist = calc_norm(this->list_nodes[i].x,this->list_nodes[i].y,this->list_nodes[i].z,\
                                p.x,p.y,p.z);

        if (dist < closest_dist)
        {
            closest_dist = dist;
            closest_index = i;
        }
    }
    return closest_index;
}

void Graph::write_info ()
{
    char filename[200];
    
    std::vector<double> segments;
    get_segments(segments);

    std::vector<double> angles;
    get_bifurcation_angles(angles);

    std::vector<double> edges;
    get_edges(edges);

    sprintf(filename,"outputs/segment_length.dat");
    write_data_to_file(filename,segments);

    sprintf(filename,"outputs/edges_length.dat");
    write_data_to_file(filename,edges);

    sprintf(filename,"outputs/bifurcation_angles.dat");
    write_data_to_file(filename,angles);

    double mean_segment_length, std_segment_length;
    compute_mean_std(segments,mean_segment_length,std_segment_length);
    printf("Segment length = %g +/- %g\n",mean_segment_length,std_segment_length);
    printf("Number of segments = %u\n",segments.size());

    double mean_edges_length, std_edges_length;
    compute_mean_std(edges,mean_edges_length,std_edges_length);
    printf("Edges length = %g +/- %g\n",mean_edges_length,std_edges_length);
    printf("Number of edges = %u\n",edges.size());

    double mean_biff_angle, std_biff_angle;
    compute_mean_std(angles,mean_biff_angle,std_biff_angle);
    printf("Bifurcation angle = %g +/- %g\n",mean_biff_angle,std_biff_angle);
    printf("Number of angles = %u\n",angles.size());
    

}
