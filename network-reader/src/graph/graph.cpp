#include "graph.h"

Graph::Graph ()
{
    total_nodes = 0;
    total_edges = 0;
}

Graph::Graph (const char filename[])
{
    // Test file format
    bool is_vtk = check_file_extension(filename,"vtk");
    bool is_txt = check_file_extension(filename,"txt");
    if (!is_vtk && !is_txt)
    {
        fprintf(stderr,"[-] ERROR! Input file must be in '.vtk' or '.txt' file format\n");
        exit(EXIT_FAILURE);
    }

    if (is_vtk)
        read_graph_from_vtk(filename);
    else if (is_txt)
        read_graph_from_txt(filename);
}

void Graph::read_graph_from_txt (const char filename[])
{
    FILE *file = fopen(filename,"r");
    if (!file)
    {
        fprintf(stderr,"[-] ERROR! Cannot open filename '%s'\n",filename);
        exit(EXIT_FAILURE);
    }    

    uint32_t num_nodes, num_edges, tmp;
    double pos[3], scalar_value;
    uint32_t dir[2];
    char str[200];
    
    // Read number of nodes
    fscanf(file,"%u",&num_nodes);

    // Initialize the list of nodes
    list_nodes.assign(num_nodes,Node());

    // Read the nodes data
    for (uint32_t i = 0; i < num_nodes; i++)
    {
        fscanf(file,"%lf %lf %lf",&pos[0],&pos[1],&pos[2]);
        
        list_nodes[i].set_node(i,pos);
    }

    // Read the edges data
    while (fscanf(file,"%u %u",&dir[0],&dir[1]) != EOF)
    {
        uint32_t src_id = dir[0];
        uint32_t dest_id = dir[1];
        double length = calc_norm(list_nodes[src_id].x,list_nodes[src_id].y,list_nodes[src_id].z,\
                                list_nodes[dest_id].x,list_nodes[dest_id].y,list_nodes[dest_id].z);

        Edge src_edge(dest_id,length);
        //Edge dest_edge(src_id,length);                        // Return edge
        list_nodes[src_id].list_edges.push_back(src_edge);
        //list_nodes[dest_id].list_edges.push_back(dest_edge);  // Return edge
    }

    fclose(file);

    // Unitary test
    assert(total_nodes == list_nodes.size());
    assert(!check_duplicates());
}

void Graph::read_graph_from_vtk (const char filename[])
{
    FILE *file = fopen(filename,"r");
    if (!file)
    {
        fprintf(stderr,"[-] ERROR! Cannot open filename '%s'\n",filename);
        exit(EXIT_FAILURE);
    }    

    uint32_t num_nodes, num_edges, tmp;
    double pos[3], scalar_value;
    uint32_t dir[2];
    char str[200];

    // Find the POINTS section
    while (fscanf(file,"%s",str) != EOF)
    {
        if (strcmp(str,"POINTS") == 0)
            break;
    }
    
    // Read number of nodes
    fscanf(file,"%u %s",&num_nodes,str);

    // Initialize the list of nodes
    list_nodes.assign(num_nodes,Node());

    // Read the nodes data
    for (uint32_t i = 0; i < num_nodes; i++)
    {
        fscanf(file,"%lf %lf %lf",&pos[0],&pos[1],&pos[2]);
        
        list_nodes[i].set_node(i,pos);

        total_nodes++;
    }

    // Find the LINES section
    while (fscanf(file,"%s",str) != EOF)
    {
        if (strcmp(str,"LINES") == 0)
            break;
    }

    // Read the edges data
    fscanf(file,"%u %u",&num_edges,&tmp);
    for (uint32_t i = 0; i < num_edges; i++)
    {
        fscanf(file,"%u %u %u",&tmp,&dir[0],&dir[1]);
        
        uint32_t src_id = dir[0];
        uint32_t dest_id = dir[1];
        double length = calc_norm(list_nodes[src_id].x,list_nodes[src_id].y,list_nodes[src_id].z,\
                                list_nodes[dest_id].x,list_nodes[dest_id].y,list_nodes[dest_id].z);

        Edge src_edge(dest_id,length);
        list_nodes[src_id].list_edges.push_back(src_edge);
        total_edges++;

        Edge dest_edge(src_id,length);                        // Return edge
        list_nodes[dest_id].list_edges.push_back(dest_edge);  // Return edge
        //total_edges++;
    }

    fclose(file);

    // Unitary test
    assert(total_nodes == list_nodes.size());
    assert(!check_duplicates());
}

void Graph::print ()
{
    printf("======================= PRINTING GRAPH ================================\n");
    for (uint32_t i = 0; i < list_nodes.size(); i++)
	{
        printf("|| %u (%g %g %g) [%u] {%d} ||",list_nodes[i].index,list_nodes[i].x,list_nodes[i].y,list_nodes[i].z,list_nodes[i].list_edges.size(),list_nodes[i].is_pmj);
        for (uint32_t j = 0; j < list_nodes[i].list_edges.size(); j++)
        {
			printf(" --> || %u %g ||",list_nodes[i].list_edges[j].dest_index,list_nodes[i].list_edges[j].length);
		}
		printf("\n");
	}
	printf("=======================================================================\n");
    printf("Number of nodes = %d\n",total_nodes);
    printf("Number of edges = %d\n",total_edges);
}

void Graph::write_network_info ()
{
    FILE *file = fopen("outputs/segments_length.dat","w+");
    for (uint32_t i = 0; i < list_nodes.size(); i++)
	{
        for (uint32_t j = 0; j < list_nodes[i].list_edges.size(); j++)
        {
			double length = list_nodes[i].list_edges[j].length;
            fprintf(file,"%g\n",length);
		}
	}
    fclose(file);

    file = fopen("outputs/bifurcations_angle.dat","w+");
    for (uint32_t i = 0; i < list_nodes.size(); i++)
	{
        if (is_bifurcation(list_nodes[i]))
        {
            uint32_t src_id = i;
            uint32_t dest_id_1 = list_nodes[i].list_edges[0].dest_index;
            uint32_t dest_id_2 = list_nodes[i].list_edges[1].dest_index;

            double u1[3], u2[3];
            build_unitary_vector(u1,src_id,dest_id_1);
            build_unitary_vector(u2,src_id,dest_id_2);

            double angle = calc_angle_between_vectors(u1,u2);
            fprintf(file,"%g\n",angle);
        }
	}
    fclose(file);
}

void Graph::set_terminals ()
{
    number_of_terminals = 0;
    for (uint32_t i = 0; i < list_nodes.size(); i++)
    {
        if (list_nodes[i].list_edges.size() == 1)
        {
            list_nodes[i].is_pmj = true;
            number_of_terminals++;
        }
    }
}

void Graph::write_terminals ()
{
    set_terminals();

    FILE *file = fopen("outputs/terminals.vtk","w+");
    fprintf(file,"# vtk DataFile Version 3.0\n");
    fprintf(file,"Terminals\n");
    fprintf(file,"ASCII\n");
    fprintf(file,"DATASET POLYDATA\n");
    fprintf(file,"POINTS %u float\n",number_of_terminals);
    for (uint32_t i = 0; i < list_nodes.size(); i++)
    {
        if (list_nodes[i].list_edges.size() == 1)
        {
            fprintf(file,"%g %g %g\n",list_nodes[i].x,list_nodes[i].y,list_nodes[i].z);
        }
    }
    fprintf(file,"VERTICES %u %u\n",number_of_terminals,number_of_terminals*2);
    for (uint32_t i = 0; i < number_of_terminals; i++)
        fprintf(file,"1 %u\n",i);

    fclose(file);
}

void Graph::build_unitary_vector (double d[], const uint32_t src_id, const uint32_t dest_id)
{
    double norm = calc_norm(list_nodes[src_id].x,list_nodes[src_id].y,list_nodes[src_id].z,\
                        list_nodes[dest_id].x,list_nodes[dest_id].y,list_nodes[dest_id].z);
    if (norm < 1.0e-08)
        norm = 1.0e-08;

    d[0] = (list_nodes[dest_id].x - list_nodes[src_id].x) / norm;
    d[1] = (list_nodes[dest_id].y - list_nodes[src_id].y) / norm;
    d[2] = (list_nodes[dest_id].z - list_nodes[src_id].z) / norm;
}

void Graph::depth_first_search (const uint32_t src_id)
{
    std::vector<bool> dfs_num;
    dfs_num.assign(total_nodes,false);

    dfs(src_id,dfs_num);

    for (uint32_t i = 0; i < dfs_num.size(); i++)
    {
        if (dfs_num[i] == DFS_WHITE)
        {
            printf("%u\n",i);
        }
    }
}

void Graph::dfs (const uint32_t src_id, std::vector<bool> &dfs_num)
{
    dfs_num[src_id] = true;
    //printf("[depth_first_search] Current on index %d\n",src_id);

    for (uint32_t i = 0; i < list_nodes[src_id].list_edges.size(); i++)
    {
        uint32_t dest_id = list_nodes[src_id].list_edges[i].dest_index;
        if (dfs_num[dest_id] == false)
        {
            dfs(dest_id,dfs_num);
        }
    }
}

bool Graph::check_duplicates ()
{
    bool ret = false;
    for (uint32_t i = 0; i < list_nodes.size(); i++)
    {
        for (uint32_t j = 0; j < list_nodes.size(); j++)
        {
            if (i != j)
            {
                double dist = calc_norm(list_nodes[i].x,list_nodes[i].y,list_nodes[i].z,\
                            list_nodes[j].x,list_nodes[j].y,list_nodes[j].z);
                if (dist < TOLERANCE_DUPLICATE)
                {
                    fprintf(stderr,"[-] ERROR! Nodes '%u' and '%u' are duplicates!\n",i,j);
                    ret = true;
                }
            }
        }
    }
    return ret;
}

bool Graph::is_bifurcation (Node u)
{
    return (u.list_edges.size() == 2) ? true : false;
}

/*
void Graph::free_list_edges (Node *node)
{
    Edge *e1 = node->list_edges;
    Edge *e2 = node->list_edges->next;

    while (e1 != NULL)
    {
        e1->next = NULL;
        e1->previous = NULL;
        e1->dest = NULL;
        delete e1;
        
        e1 = e2;
        if (e2 != NULL)
            e2 = e2->next;
    }

    node->list_edges = NULL;
    node->last_edge = NULL;
    node->num_edges = 0;
}

void Graph::free_list_nodes ()
{
    Node *n1 = list_nodes;
    Node *n2 = list_nodes->next;

    while (n1 != NULL)
    {
        if (n1->list_edges)
            free_list_edges(n1);
        
        n1->next = NULL;
        n1->previous = NULL;
        delete n1;

        n1 = n2;
        if (n2 != NULL)
            n2 = n2->next;

    }

    list_nodes = NULL;
    last_node = NULL;
    total_nodes = 0;
    total_edges = 0;
}

Graph::~Graph ()
{
    if (list_nodes)
        free_list_nodes();
}

void Graph::insert_node_graph (const double pos[])
{
    // TODO: Add a check for duplicates ...

    Node *node = new Node(total_nodes,pos);

    // First node of the list
    if (last_node == NULL)
    {
        list_nodes = node;
        last_node = node;
    }
    // Insert after the last node and update teh pointer
    else
    {
        last_node->next = node;
        node->previous = last_node;
        last_node = node;
    }
    total_nodes++;

}

void Graph::insert_edge_graph (const int id_1, const int id_2)
{
    // TODO: Check if the edge is invalid
    //if (id_1 == id_2) return;

    Node *n1 = search_node(id_1);
    Node *n2 = search_node(id_2);

    double norm = calc_norm(n1->x,n1->y,n1->z,n2->x,n2->y,n2->z);
    Edge *edge = new Edge(id_2,norm,n2);

    // First edge
    if (n1->last_edge == NULL)
    {
        n1->list_edges = edge;
        n1->last_edge = edge;
    }
    // Insert after the last edge from the list 
    else
    {
        n1->last_edge->next = edge;
        edge->previous = n1->last_edge;
        n1->last_edge = edge;
    }
    
    // Increment the number of edges from the origin Node
    n1->num_edges++;
    
    // Increment the total number of edges of the graph
    total_edges++;
}

// TODO: Improve this code ...
void Graph::remove_node_graph (const int id)
{
    // 1) Remove any other edge that points to the deleted node
    Node *node_to_delete = NULL;
    Node *ptr = list_nodes;

    while (ptr != NULL)
    {
        // Save the pointer of the node to be deleted
        if (ptr->index == id)
            node_to_delete = ptr;
        
        remove_edge_graph(ptr->index,id);
        ptr = ptr->next;
    }

    // 2) Remove all the edges from the node to be deleted
    while (node_to_delete->num_edges != 0)
        remove_edge_graph(node_to_delete->index,node_to_delete->list_edges->index);
    
    // 3) Set the pointers based on each case
    // Case 1: First node
    if (node_to_delete == list_nodes)
    {
        list_nodes = node_to_delete->next;
        list_nodes->previous = NULL;
        node_to_delete->next = NULL;
    }
    // Case 2: Last node
    else if (node_to_delete == last_node)
    {
        last_node = node_to_delete->previous;
        last_node->next = NULL;
        node_to_delete->previous = NULL;
    }
    // Case 3: Middle node
    else
    {
        node_to_delete->previous->next = node_to_delete->next;
        node_to_delete->next->previous = node_to_delete->previous;
        node_to_delete->next = NULL;
        node_to_delete->previous = NULL;
    }
    delete node_to_delete;
    total_nodes--;

    // Corner case
    if (total_nodes == 0)
        last_node = NULL;

    // 4) Update the indexes of the graph
    ptr = list_nodes;
    while (ptr != NULL)
    {
        Edge *ptrl = ptr->list_edges;

        // The nodes greater than the deleted one need to be updated ...
        if (ptr->index > id)
            ptr->index--;
        
        while (ptrl != NULL)
        {
            // To the same for teh index of the edges ...
            if (ptrl->index > id)
                ptrl->index--;
            ptrl = ptrl->next;
        }

        ptr = ptr->next;
    }

}

void Graph::remove_edge_graph (const int id_1, const int id_2)
{
    Node *node = search_node(id_1);
    
    Edge *ptrl = node->list_edges;
    while (ptrl != NULL)
    {
        if (ptrl->index == id_2)
        {
            // Remove edge
            // Case 1: First edge
            if (ptrl == node->list_edges)
            {
                node->list_edges = ptrl->next;
                ptrl->previous = NULL;
                ptrl->next = NULL;
                ptrl->dest = NULL;

            }
            // Case 2: Last edge
            else if (ptrl == node->last_edge)
            {
                node->last_edge = ptrl->previous;
                node->last_edge->next = NULL;
                ptrl->previous = NULL;
                ptrl->dest = NULL;

            }
            // Case 3: Middle edge
            else
            {
                ptrl->previous->next = ptrl->next;
                ptrl->next->previous = ptrl->previous;
                ptrl->next = NULL;
                ptrl->previous = NULL;

            }
            delete ptrl;
            node->num_edges--;
            total_edges--;

            // Corner case
            if (node->num_edges == 0)
                node->last_edge = NULL;

            return;
        }
        ptrl = ptrl->next;
    }

    //fprintf(stderr,"[-] ERROR! Edge (%d,%d) not found!\n",id_1,id_2);
    //exit(EXIT_FAILURE);
}

Node* Graph::search_node (const int id)
{
    Node *tmp = list_nodes;
    while (tmp != NULL)
    {
        if (tmp->index == id)
            return tmp;
        tmp = tmp->next;
    }
    //fprintf(stderr,"[-] ERROR! Node %d not found!\n",id);

    return NULL;
}

Node::Node (const int id, const double pos[])
{
	index = id;

	x = pos[0];
	y = pos[1];
	z = pos[2];

	next = NULL;
    previous = NULL;

    num_edges = 0;
	list_edges = NULL;
    last_edge = NULL;

    is_pmj = false;
}

Edge::Edge (const int id, const double l, Node *destination)
{
    index = id;
    length = l;
    dest = destination;

    next = NULL;
    previous = NULL;
}

void Graph::print ()
{
    Node *ptr = list_nodes;
	printf("======================= PRINTING GRAPH ================================\n");
	while (ptr != NULL)
	{
        printf("|| %d (%.2lf %.2lf %.2lf) Edges = %d ||",ptr->index,ptr->x,ptr->y,ptr->z,ptr->num_edges);
        Edge *ptrl = ptr->list_edges;
		while (ptrl != NULL)
		{
			printf(" --> || %d %.2lf ||",ptrl->index,ptrl->length);
			ptrl = ptrl->next;
		}
		printf("\n");
		ptr = ptr->next;
	}
	printf("=======================================================================\n");
    printf("Number of nodes = %d\n",total_nodes);
    printf("Number of edges = %d\n",total_edges);
}

void Graph::write_VTK (const char filename[])
{
    Node *ptr;
    Edge *ptrl;

    FILE *file = fopen(filename,"w+");

    fprintf(file,"# vtk DataFile Version 3.0\n");
    fprintf(file,"Graph\n");
    fprintf(file,"ASCII\n");
    fprintf(file,"DATASET POLYDATA\n");
    fprintf(file,"POINTS %d float\n",total_nodes);
    ptr = list_nodes;
    while (ptr != NULL)
    {
        fprintf(file,"%g %g %g\n",ptr->x,ptr->y,ptr->z);
        ptr = ptr->next;
    }
    fprintf(file,"LINES %d %d\n",total_edges,total_edges*3);
    ptr = list_nodes;
    while (ptr != NULL)
    {
        ptrl = ptr->list_edges;
        while (ptrl != NULL)
        {
            fprintf(file,"2 %d %d\n",ptr->index,ptrl->index);
            ptrl = ptrl->next;
        }
        ptr = ptr->next;
    }

    fclose(file);
}

void Graph::write_pmj_config_file (const char filename[])
{
    const double RPMJ = 1.0e+03;

    FILE *file = fopen(filename,"w+");

    Node *tmp = this->list_nodes;
    while (tmp != NULL)
    {
        if (is_terminal(tmp))
        {
            if (tmp->is_pmj)
                fprintf(file,"%u %g\n",tmp->index,RPMJ);
            else
                fprintf(file,"%u %g\n",tmp->index,__DBL_MAX__);
        }
            
        tmp = tmp->next;
    }

    fclose(file);
}

bool is_terminal (Node *u)
{
    if (u->num_edges == 0 && u->index != 0)
        return true;
    else
        return false;
}

bool is_bifurcation (Node *u)
{
    if (u->num_edges == 2)
        return true;
    else
        return false;
}

void Graph::depth_first_search ()
{
    uint32_t source_index = 142;
    Node *source_node = search_node(source_index);

    vector<int> dfs_num;
    dfs_num.assign(total_nodes,DFS_WHITE);

    dfs(source_node,dfs_num);

    for (uint32_t i = 0; i < dfs_num.size(); i++)
    {
        if (dfs_num[i] == DFS_WHITE)
        {
            printf("%u\n",i);
        }
    }
}

void Graph::dfs (Node *u, vector<int> &dfs_num)
{
    int u_index = u->index;
    Edge *v = u->list_edges;

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

void Graph::breadth_first_search ()
{
    int *parents = new int[total_nodes]();

    int source_index = 0;
    //int source_index = 142;
    Node *source_node = search_node(source_index);
    parents[source_index] = -1;

    map<int,int> dist;              // Distance from source to the other nodes
    dist[source_index] = 0;         // Distance from source to source is zero

    queue<Node*> q;
    q.push(source_node);            // Enqueue the source first

    while (!q.empty())
    {
        Node *u = q.front(); q.pop();
        int u_index = u->index;
        //printf("[breadth_first_search] Current on index %d\n",u_index);

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

    int max_level = -1;
    for (auto it = dist.begin(); it != dist.end(); ++it)
    {
        //printf("[breadth_first_search] Node %d -- Distance from source %d\n",it->first,it->second);

        if (it->second > max_level)
            max_level = it->second;
    }

    vector<int> nodes_per_level;
    nodes_per_level.assign(max_level+1,0);
    for (auto it = dist.begin(); it != dist.end(); ++it)
    {
        int node_index = it->first;
        int level = it->second;

        nodes_per_level[level]++;
    }

    for (uint32_t i = 0; i < nodes_per_level.size(); i++)
    {
        printf("[breadth_first_search] Level %3d ",i);
        print_stars(nodes_per_level[i]);
        printf(" (%d)\n",nodes_per_level[i]);
    }

    // Discover the branch with the longest size
    double max_size = __DBL_MIN__;
    int max_num_biff = 0;
    int max_index = -1;
    Node *tmp = this->list_nodes;
    while (tmp != NULL)
    {
        double size = 0.0;
        int num_biff = 0;
        if (is_terminal(tmp))
        {
            size = calculate_branch_size(parents,tmp->index,num_biff);
            if (num_biff > max_num_biff)
            {
                max_size = size;
                max_index = tmp->index;
                max_num_biff = num_biff;
            }
        }
             
        tmp = tmp->next;
    }

    write_longest_segment(parents,max_index);
    printf("Longest segment %d -- Number of bifurcations %d -- Size = %g\n",max_index,max_num_biff,max_size);

    delete [] parents;
}

double Graph::calculate_branch_size (const int parents[], const int ref_index, int &num_biff)
{
    double size = 0.0;
    int counter = 0;
    int current_index = ref_index;
    int prev_index = ref_index;
    while (current_index != -1)
    {

        prev_index = current_index;
        current_index = parents[current_index];

        Node *prev = search_node(prev_index);
        Node *cur = search_node(current_index);

        if (prev && cur)
            size += calc_norm(prev->x,prev->y,prev->z,cur->x,cur->y,cur->z);
        
        if (prev->num_edges == 2)
            counter++;
        
    }
    
    num_biff = counter;

    return size;
}

void Graph::write_longest_segment (const int parents[], const int ref_index)
{
    int counter = 0;
    int current_index = ref_index;
    int prev_index = ref_index;
    while (current_index != -1)
    {
        //printf("Node %d -- Parent %d\n",current_index,parents[current_index]);
        current_index = parents[current_index];
        counter++;
    }

    char filename[200];
    sprintf(filename,"outputs/longest_segment_%d.vtk",ref_index);
    FILE *file = fopen(filename,"w+");

    fprintf(file,"# vtk DataFile Version 3.0\n");
    fprintf(file,"Graph\n");
    fprintf(file,"ASCII\n");
    fprintf(file,"DATASET POLYDATA\n");
    fprintf(file,"POINTS %d float\n",counter);
    
    
    current_index = ref_index;
    while (current_index != -1)
    {
        prev_index = current_index;
        current_index = parents[current_index];

        Node *prev = search_node(prev_index);
        fprintf(file,"%g %g %g\n",prev->x,prev->y,prev->z);    

    }

    fprintf(file,"LINES %d %d\n",counter-1,(counter-1)*3);    
    for (uint32_t i = 0; i < (counter-1); i++)
        fprintf(file,"2 %d %d\n",i,i+1);

    fclose(file);

}

void Graph::check_duplicates ()
{

    // Check duplicates
    Node *tmp = this->get_list_nodes();
    while (tmp != NULL)
    {
        Node *tmp2 = this->get_list_nodes();
        while (tmp2 != NULL)
        {
            if (tmp->x == tmp2->x && tmp->y == tmp2->y && tmp->z == tmp2->z && tmp->index != tmp2->index)
            {
                printf("[purkinje] Duplicates are indexes: %u and %u --> (%g,%g,%g) x (%g,%g,%g)\n",tmp->index,tmp2->index,tmp->x,tmp->y,tmp->z,tmp2->x,tmp2->y,tmp2->z);
                printf("\t|| %u has %u edges || %u has %u edges ||\n",tmp->index,tmp->num_edges,tmp2->index,tmp2->num_edges);
                remove_node_graph(tmp2->index);
            }
            tmp2 = tmp2->next;
        }
        tmp = tmp->next;
    }
}

void Graph::build_unitary_vector (Node *u, Node *v, double d[])
{
    double norm = calc_norm(u->x,u->y,u->z,v->x,v->y,v->z);
    if (norm < 1.0e-08)
        norm = 1.0e-08;

    d[0] = (v->x - u->x) / norm;
    d[1] = (v->y - u->y) / norm;
    d[2] = (v->z - u->z) / norm;
}

double Graph::calc_angle_between_vectors (const double u[], const double v[])
{
    double dot_product = u[0]*v[0] + u[1]*v[1] + u[2]*v[2];

    double angle_radians = acos(dot_product);

    // Return the angle in degrees
    return angle_radians * 180.0 / M_PI;
}

void Graph::write_network_info ()
{
    // [SEGMENTS]
    FILE *file = fopen("outputs/segments_length.dat","w+");
    Node *u = this->get_list_nodes();
    while (u != NULL)
    {
        Edge *v = u->list_edges;
        while (v != NULL)
        {
            fprintf(file,"%g\n",v->length);
            v = v->next;
        }

        u = u->next;
    }
    fclose(file);

    file = fopen("outputs/bifurcations_angle.dat","w+");
    u = this->get_list_nodes();
    while (u != NULL)
    {
        if (is_bifurcation(u))
        {
            Edge *v1 = u->list_edges;
            Edge *v2 = u->list_edges->next;

            double u1[3], u2[3];
            build_unitary_vector(u,v1->dest,u1);
            build_unitary_vector(u,v2->dest,u2);

            double angle = calc_angle_between_vectors(u1,u2);
            fprintf(file,"%g\n",angle);
        }
        u = u->next;
    }
    fclose(file);
}

void Graph::print_terminals ()
{
    uint32_t counter = 0;
    Node *tmp = this->get_list_nodes();
    while (tmp != NULL)
    {
        if (is_terminal(tmp))
        {
            printf("%g %g %g\n",tmp->x,tmp->y,tmp->z);
            counter++;
        }
            
        tmp = tmp->next;
    }
    printf("%u\n",counter);
}

void print_stars (const int number)
{
    for (int i = 0; i < number; i++)
        printf("*"); 
}


bool Graph::is_duplicate (const double pos[])
{
	Node *ptr = list_nodes;
	while (ptr != NULL)
	{	
		if (calc_norm(pos[0],pos[0],pos[0],ptr->x,ptr->y,ptr->z) < TOLERANCE_DUPLICATE)
			return true;
		ptr = ptr->next;
	}
	return false;
}

void Graph::dijkstra (int s)
{
    printf("[!] Running Dijkstra ... \n");

    if (!dist)
    {
        dist = (double*)malloc(sizeof(double)*total_nodes);
    }
        
    for (int i = 0; i < total_nodes; i++) 
        dist[i] = INF;
    dist[s] = 0;

    priority_queue< pair<double,int>, vector< pair<double,int> >, greater< pair<double,int> > > pq;
    pq.push(make_pair(0,s));

    while (!pq.empty())
    {
        pair<double,int> front = pq.top(); pq.pop();
        double d = front.first;
        int u = front.second;
        if (d > dist[u]) 
            continue;
        Edge *ptrl = search_node(u)->list_edges;
        while (ptrl != NULL)
        {
            int id = ptrl->id;
            double w = ptrl->w; 
            if (dist[u] + w < dist[id])
            {
                dist[id] = dist[u] + w;
                pq.push(make_pair(dist[id],id));
            }
            ptrl = ptrl->next;
        }
    }

}

*/
