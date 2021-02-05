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

    // Run a Dijkstra shortest path to fill the 'dist' and 'parent' arrays
    dijkstra(0);

    // Compute the LAT of all cells
    compute_activation_times(REF_CV);

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

uint32_t Graph::get_closest_terminal_point (Point p, const bool is_reference)
{
    uint32_t biff_value = (is_reference) ? 1 : 0;
    uint32_t closest_index = 0;
    double closest_dist = __DBL_MAX__;
    for (uint32_t i = 0; i < this->list_nodes.size(); i++)
    {
        if (this->list_nodes[i].list_edges.size() == biff_value)
        {
            double dist = calc_norm(this->list_nodes[i].x,this->list_nodes[i].y,this->list_nodes[i].z,\
                                p.x,p.y,p.z);

            if (dist < closest_dist)
            {
                closest_dist = dist;
                closest_index = i;
            }
        }
    }
    return closest_index;
}

// Using C++/STD library
void Graph::dijkstra (const uint32_t src_id)
{
    uint32_t np = this->total_nodes;
    dist.assign(np,__DBL_MAX__);
    parent.assign(np,-1);
    dist[src_id] = 0.0;

    std::priority_queue< std::pair<double,uint32_t>, std::vector< std::pair<double,uint32_t> >, std::greater< std::pair<double,uint32_t> > > pq;
    pq.push(std::make_pair(0.0,src_id));

    while (!pq.empty())
    {
        std::pair<double,uint32_t> front = pq.top(); pq.pop();
        double d = front.first;
        uint32_t u = front.second;
        if (d > dist[u]) 
            continue;
        
        for (uint32_t i = 0; i < this->list_nodes[u].list_edges.size(); i++)
        {
            uint32_t v = this->list_nodes[u].list_edges[i].dest_id;
            double w = this->list_nodes[u].list_edges[i].length;

            if (dist[u] + w < dist[v])
            {
                dist[v] = dist[u] + w;
                parent[v] = u;
                pq.push(std::make_pair(dist[v],v));
            }
        }
    }

    //for (uint32_t i = 0; i < dist.size(); i++)
    //    printf("Node %u -- Parent = %d -- Dist = %g\n",i,parent[i],dist[i]);
}

// Using PQUEUE library
void Graph::dijkstra_2 (const uint32_t src_id)
{
    // Initialize the shortest distance array
    uint32_t np = this->total_nodes;
    dist.assign(np,__DBL_MAX__);
    dist[src_id] = 0.0;

    pqueue_t *pq;
	node_t   *ns;
	node_t   *n;

    // Initialize the priority queue
	ns = (struct node_t*)malloc(np * sizeof(node_t));
	pq = pqueue_init(np, cmp_pri, get_pri, set_pri, get_pos, set_pos);
	if (!(ns && pq))
    {
        fprintf(stderr,"ERROR! Bad alloc!\n");
        exit(EXIT_FAILURE); 
    }
    ns[src_id].pri = 0.0; ns[src_id].val = src_id; pqueue_insert(pq, &ns[src_id]);

    while ((n = (node_t*)pqueue_pop(pq)))
    {
        double d = n->pri;
        uint32_t u = n->val;
        if (d > dist[u]) 
            continue;
        
        for (uint32_t i = 0; i < this->list_nodes[u].list_edges.size(); i++)
        {
            uint32_t v = this->list_nodes[u].list_edges[i].dest_id;
            double w = this->list_nodes[u].list_edges[i].length;

            if (dist[u] + w < dist[v])
            {
                dist[v] = dist[u] + w;
                ns[v].pri = dist[v]; ns[v].val = v; pqueue_insert(pq, &ns[v]);
            }
        }
    }

	pqueue_free(pq);
	free(ns);

    //for (uint32_t i = 0; i < dist.size(); i++)
    //    printf("Node %u -- Dist = %g\n",i,dist[i]);
}

// Purkinje propagation velocity {um/ms}
void Graph::calculate_activation_time (std::vector<Point> pmj_points, const double cv)
{
    
}

void Graph::fill_terminal_indexes (std::vector<Point> pmj_points, const bool is_reference)
{
    for (uint32_t i = 0; i < pmj_points.size(); i++)
    {
        uint32_t index = get_closest_terminal_point(pmj_points[i],is_reference);
        this->terminals_indexes.push_back(index);
        //printf("%u %g %g %g\n",index,this->list_nodes[index].x,this->list_nodes[index].y,this->list_nodes[index].z);
    }
    char filename[200];
    if (is_reference)
        sprintf(filename,"outputs/reference_terminals.vtk");
    else
        sprintf(filename,"outputs/aproximation_terminals.vtk");
    //write_terminals(filename);
}

void Graph::compute_errors (Graph *input, std::string pmj_filename)
{
    std::vector<Point> pmj_points;
    read_points(pmj_filename,pmj_points);

    this->fill_terminal_indexes(pmj_points,true);
    input->fill_terminal_indexes(pmj_points,false);

/*
    uint32_t ref_min_lat_id, ref_max_lat_id;
    double ref_min_lat, ref_max_lat;
    uint32_t aprox_min_lat_id, aprox_max_lat_id;
    double aprox_min_lat, aprox_max_lat;
    double percentage_2ms, percentage_5ms;
    double rmse, rrmse, max_error;

    std::vector<double> ref_lat = this->lat;
    std::vector<double> aprox_lat = input->lat;
    std::vector<uint32_t> ref_term_ids = this->terminals_indexes;
    std::vector<uint32_t> aprox_term_ids = input->terminals_indexes;
    
    compute_rmse_rrmse_maxerror(ref_lat,aprox_lat,ref_term_ids,aprox_term_ids,rmse,rrmse,max_error);
    compute_min_max_lat(ref_lat,ref_term_ids,ref_min_lat,ref_max_lat,ref_min_lat_id,ref_max_lat_id);
    compute_min_max_lat(aprox_lat,aprox_term_ids,aprox_min_lat,aprox_max_lat,aprox_min_lat_id,aprox_max_lat_id);
    compute_epsilon_percentage(ref_lat,aprox_lat,ref_term_ids,aprox_term_ids,2,percentage_2ms);
    compute_epsilon_percentage(ref_lat,aprox_lat,ref_term_ids,aprox_term_ids,5,percentage_5ms);

    double ref_dist_min_lat = this->dist[ref_min_lat_id];
    double ref_dist_max_lat = this->dist[ref_max_lat_id];
    double aprox_dist_min_lat = input->dist[aprox_min_lat_id];
    double aprox_dist_max_lat = input->dist[aprox_max_lat_id];
*/
    //printf("========================================================================================================================================================================\n");
    /*
    printf("Max error PMJ's = %g ms\n",max_error);
    printf("RMSE PMJ's = %g ms\n",rmse);
    printf("RRMSE PMJ's = %g %%\n",rrmse*100.0);
    printf("Epsilon < 2ms = %g %%\n",percentage_2ms);
    printf("Epsilon < 5ms = %g %%\n",percentage_5ms);
    printf("[reference] Min. LAT PMJ's = %g || Max. LAT PMJ's = %g || Min. LAT cell id = %u || Max. LAT cell id = %u || Min. LAT dist = %g || Max. LAT dist = %g\n",ref_min_lat,ref_max_lat,ref_min_lat_id,ref_max_lat_id,ref_dist_min_lat,ref_dist_max_lat);
    printf("[aproximation] Min. LAT PMJ's = %g || Max. LAT PMJ's = %g || Min. LAT cell id = %u || Max. LAT cell id = %u || Min. LAT dist = %g || Max. LAT dist = %g\n",aprox_min_lat,aprox_max_lat,aprox_min_lat_id,aprox_max_lat_id,aprox_dist_min_lat,aprox_dist_max_lat);
    */
    //printf("%.2lf\t%.2lf\t%.2lf\t%.2lf\t%.2lf\t%.2lf\t%.2lf\n",aprox_min_lat,aprox_max_lat,max_error,rmse,rrmse*100.0,percentage_2ms,percentage_5ms);
    //printf("========================================================================================================================================================================\n");

    // [GLOBAL] Recalculate the propagation velocity of the aproximation Purkinje network using proportion calculus
    //double new_cv = aprox_dist_max_lat / ref_max_lat;
    //double new_diameter = calculate_diameter(new_cv);
    //input->compute_activation_times(new_cv);
    //printf("New CV = %g um/mm\n",new_cv);
    //printf("New diameter = %g um\n",new_diameter);
    //printf("\n");
/*
    std::vector<double> aprox_lat_2 = input->lat;
    compute_rmse_rrmse_maxerror(ref_lat,aprox_lat_2,ref_term_ids,aprox_term_ids,rmse,rrmse,max_error);
    compute_min_max_lat(aprox_lat_2,aprox_term_ids,aprox_min_lat,aprox_max_lat,aprox_min_lat_id,aprox_max_lat_id);
    compute_epsilon_percentage(ref_lat,aprox_lat_2,ref_term_ids,aprox_term_ids,2,percentage_2ms);
    compute_epsilon_percentage(ref_lat,aprox_lat_2,ref_term_ids,aprox_term_ids,5,percentage_5ms);
*/
    //printf("========================================================================================================================================================================\n");
    /*
    printf("Max error PMJ's = %g ms\n",max_error);
    printf("RMSE PMJ's = %g ms\n",rmse);
    printf("RRMSE PMJ's = %g %%\n",rrmse*100.0);
    printf("Epsilon < 2ms = %g %%\n",percentage_2ms);
    printf("Epsilon < 5ms = %g %%\n",percentage_5ms);
    printf("[reference] Min. LAT PMJ's = %g || Max. LAT PMJ's = %g || Min. LAT cell id = %u || Max. LAT cell id = %u || Min. LAT dist = %g || Max. LAT dist = %g\n",ref_min_lat,ref_max_lat,ref_min_lat_id,ref_max_lat_id,ref_dist_min_lat,ref_dist_max_lat);
    printf("[aproximation] Min. LAT PMJ's = %g || Max. LAT PMJ's = %g || Min. LAT cell id = %u || Max. LAT cell id = %u || Min. LAT dist = %g || Max. LAT dist = %g\n",aprox_min_lat,aprox_max_lat,aprox_min_lat_id,aprox_max_lat_id,aprox_dist_min_lat,aprox_dist_max_lat);
    */
    //printf("%.2lf\t%.2lf\t%.2lf\t%.2lf\t%.2lf\t%.2lf\t%.2lf\n",aprox_min_lat,aprox_max_lat,max_error,rmse,rrmse*100.0,percentage_2ms,percentage_5ms);
    //printf("========================================================================================================================================================================\n");

}

void Graph::compute_epsilon_percentage (std::vector<double> ref_lat, std::vector<double> aprox_lat,\
                                    std::vector<uint32_t> ref_term_ids, std::vector<uint32_t> aprox_term_ids,\
                                    const double epsilon, double &percentage)
{
    uint32_t counter = 0;
    uint32_t n = ref_term_ids.size();
    for (uint32_t i = 0; i < n; i++)
    {
        uint32_t ref_index = ref_term_ids[i];
        uint32_t aprox_index = aprox_term_ids[i];
        double error = fabs(ref_lat[ref_index] - aprox_lat[aprox_index]);
        if (error < epsilon)
            counter++;
    }    
    percentage = (double)counter / (double)n * 100.0;
}

void Graph::compute_activation_times (const double cv)
{
    this->lat.clear();
    for (uint32_t i = 0; i < this->total_nodes; i++)
    {
        this->lat.push_back(this->dist[i]/cv);
    }
}

void Graph::compute_min_max_lat (std::vector<double> lat, std::vector<uint32_t> term_ids, double &min_value, double &max_value, uint32_t &min_id, uint32_t &max_id)
{
    max_value = __DBL_MIN__;
    min_value = __DBL_MAX__;
    for (uint32_t i = 0; i < term_ids.size(); i++)
    {
        uint32_t index = term_ids[i];
        if (lat[index] > max_value)
        {
            max_value = lat[index];
            max_id = index;
        }
        if (lat[index] < min_value)
        {
            min_value = lat[index];
            min_id = index;
        }
    }
}

void Graph::compute_rmse_rrmse_maxerror (std::vector<double> ref_lat, std::vector<double> aprox_lat,\
                                        std::vector<uint32_t> ref_term_ids, std::vector<uint32_t> aprox_term_ids,\
                                        double &rmse, double &rrmse, double &max_error)
{
    max_error = __DBL_MIN__;
    uint32_t n = ref_term_ids.size();
    double sum_num = 0.0;
    double sum_den = 0.0;
    for (uint32_t i = 0; i < n; i++)
    {
        uint32_t ref_index = ref_term_ids[i];
        uint32_t aprox_index = aprox_term_ids[i];
        double error = fabs(ref_lat[ref_index] - aprox_lat[aprox_index]);
        //printf("%u %g\n",i,error);
        if (error > max_error)
            max_error = error;

        sum_num += powf(error,2);
        sum_den += powf(ref_lat[ref_index],2);
    }    
    double l2_norm = sqrt(sum_den);
    rmse = sqrt(sum_num/(double)n);
    rrmse = sqrt(sum_num/sum_den);
    //printf("\n");
}

double calculate_diameter (const double cv)
{
    static const double Gi = 7.9;
    static const double Cf = 3.4;
    static const double tauf = 0.1;
    double cv_m_per_s = cv / 1000.0;

    return (cv_m_per_s*cv_m_per_s*4.0*Cf*tauf) / (Gi) * 100.0;
}

double calculate_proportion (const double ref_diameter, const double ref_max_lat, const double aprox_max_lat)
{
    return ref_diameter * aprox_max_lat / ref_max_lat;
}

double calculate_propagation_velocity (const double diameter)
{
    static const double Gi = 7.9;
    static const double Cf = 3.4;
    static const double tauf = 0.1;

    return sqrt( (Gi*diameter)/(4.0*Cf*tauf) ) * 0.1;
}

double adjust_propagation_velocity (const double ref_cv, const double ref_max_lat, const double aprox_max_lat, double &new_diameter)
{
    double ref_diameter = calculate_diameter(ref_cv);
    new_diameter = calculate_proportion(ref_diameter,ref_max_lat,aprox_max_lat);
    return calculate_propagation_velocity(new_diameter)*1000.0;
}

void Graph::write_terminals (const char filename[])
{
    FILE *file = fopen(filename,"w+");
    /*
    fprintf(file,"# vtk DataFile Version 4.1\n");
    fprintf(file,"vtk output\n");
    fprintf(file,"ASCII\n");
    fprintf(file,"DATASET POLYDATA\n");
    
    fprintf(file,"POINTS %u float\n",this->terminals_indexes.size());
    for (uint32_t i = 0; i < this->terminals_indexes.size(); i++)
    {
        uint32_t index = this->terminals_indexes[i];
        fprintf(file,"%g %g %g\n",this->list_nodes[index].x,this->list_nodes[index].y,this->list_nodes[index].z);
    }
    
    fprintf(file,"VERTICES %u %u\n",this->terminals_indexes.size(),this->terminals_indexes.size()*2);
    for (uint32_t i = 0; i < this->terminals_indexes.size(); i++)
        fprintf(file,"1 %u\n",i);
    
    fprintf(file,"POINT_DATA %u\n",this->terminals_indexes.size());
    fprintf(file,"SCALARS LAT float\n");
    fprintf(file,"LOOKUP_TABLE default\n");
    for (uint32_t i = 0; i < this->terminals_indexes.size(); i++)
    {
        uint32_t index = this->terminals_indexes[i];
        fprintf(file,"%g\n",this->lat[index]);
    }
    */

    for (uint32_t i = 0; i < this->terminals_indexes.size(); i++)
    {
        uint32_t index = this->terminals_indexes[i];
        fprintf(file,"%g %g %g %g\n",this->list_nodes[index].x,this->list_nodes[index].y,this->list_nodes[index].z,this->lat[index]);
    }

    fclose(file);
}

void Graph::write_LAT (const char filename[])
{
    FILE *file = fopen(filename,"w+");
    fprintf(file,"# vtk DataFile Version 4.1\n");
    fprintf(file,"vtk output\n");
    fprintf(file,"ASCII\n");
    fprintf(file,"DATASET POLYDATA\n");
    
    fprintf(file,"POINTS %u float\n",this->list_nodes.size());
    for (uint32_t i = 0; i < this->list_nodes.size(); i++)
    {
        fprintf(file,"%g %g %g\n",this->list_nodes[i].x,this->list_nodes[i].y,this->list_nodes[i].z);
    }
    fprintf(file,"LINES %u %u\n",this->total_edges,this->total_edges*3);
    for (uint32_t i = 0; i < this->list_nodes.size(); i++)
    {
        uint32_t u = this->list_nodes[i].id;
        for (uint32_t j = 0; j < this->list_nodes[i].list_edges.size(); j++)
        {
            uint32_t v = this->list_nodes[i].list_edges[j].dest_id;
            fprintf(file,"2 %u %u\n",u,v);
        }
    }
    fprintf(file,"POINT_DATA %u\n",this->list_nodes.size());
    fprintf(file,"SCALARS LAT float\n");
    fprintf(file,"LOOKUP_TABLE default\n");
    for (uint32_t i = 0; i < this->list_nodes.size(); i++)
    {
        fprintf(file,"%g\n",this->lat[i]);
    }
    fclose(file);

    for (uint32_t i = 0; i < this->terminals_indexes.size(); i++)
        printf("Active PMJ = %u || LAT = %g\n",this->terminals_indexes[i],this->lat[this->terminals_indexes[i]]);
}
