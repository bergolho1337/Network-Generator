#include "graph.h"

Graph::Graph ()
{
    last_node = NULL;
    list_nodes = NULL;
    total_nodes = 0;
    total_edges = 0;
    dist = NULL;
}

void Graph::free_list_edges (Node *node)
{
    Edge *e1 = node->list_edges;
    Edge *e2 = node->list_edges->next;

    while (e1 != NULL)
    {
        e1->next = NULL;
        delete e1;
        e1 = e2;
        if (e2 != NULL)
            e2 = e2->next;
    }

    node->list_edges = NULL;
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
        delete n1;
        n1 = n2;
        if (n2 != NULL)
            n2 = n2->next;

    }
}

Graph::~Graph ()
{
    if (list_nodes)
        free_list_nodes();
    
    if (dist)
        free(dist);
}

Node::Node (const int i, const double pos[], const double d[])
{
	id = i;

	x = pos[0];
	y = pos[1];
	z = pos[2];

	d_ori[0] = d[0];
	d_ori[1] = d[1];
	d_ori[2] = d[2];

	is_terminal = false;

	num_edges = 0;
	next = NULL;
	list_edges = NULL;
}

Edge::Edge (const int i, const double weight, Node *destination)
{
    id = i;
    w = weight;
    dest = destination;
    next = NULL;
}

void Graph::calc_original_growth_direction (double d_ori[], const Node *prev,\
					    const double x, const double y, const double z)
{
	if (prev == NULL)
	{
		d_ori[0] = 0.0;
		d_ori[1] = 0.0;
		d_ori[2] = 0.0;
	}
	else
	{
		double norm = calc_norm(prev->x,prev->y,prev->z,x,y,z);
		d_ori[0] = (x - prev->x) / norm;
		d_ori[1] = (y - prev->y) / norm;
		d_ori[2] = (z - prev->z) / norm;
	}
}

// Insert a node in the graph by receiving its position (x,y,z) and a reference to its predecessor
Node* Graph::insert_node_graph (const double pos[], const Node *prev)
{
	// First check if the node already exists in the network
	if (!is_duplicate(pos))
	{
		double d_ori[3];
		calc_original_growth_direction(d_ori,prev,pos[0],pos[1],pos[2]);

		Node *tmp = list_nodes;		
		Node *node = new Node(total_nodes++,pos,d_ori);
		
		// First node of the list
		if (!tmp)
		{
			list_nodes = node;
			last_node = node;
		}
		// Insert after the last node and update the pointer
		else
		{
			last_node->next = node;
			last_node = last_node->next;
		}
	
		// Return a reference to the node that has been inserted
		return node;
	}
	else
	{
		printf("[-] ERROR! Position (%.10lf,%10lf,%10lf) has already been taken by a Node!\n",pos[0],pos[1],pos[2]);
		return NULL;
	}
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

Node* Graph::search_node (const int id)
{
    Node *tmp = list_nodes;
    while (tmp != NULL)
    {
        if (tmp->id == id)
            return tmp;
        tmp = tmp->next;
    }
    cerr << "[-] ERROR! Node " << id << " not found!" << endl;
    return NULL;
}

void Graph::insert_edge_graph (const int id_1, const int id_2)
{
    Node *n1, *n2;
    Edge *edge;
    double norm;

    // Check if the edge is invalid
    if (id_1 == id_2) return;

    n1 = search_node(id_1);
    n2 = search_node(id_2);

    //printf("Insert edge (%d,%d) -- (%.10lf,%.10lf,%.10lf) -> (%.10lf,%.10lf,%.10lf)\n",n1->id,n2->id,\
								n1->x,n1->y,n1->z,n2->x,n2->y,n2->z);

    norm = calc_norm(n1->x,n1->y,n1->z,n2->x,n2->y,n2->z);
    edge = new Edge(id_2,norm,n2);

    // First edge
    if (!n1->list_edges)
        n1->list_edges = edge;
    // Iterate over the list and insert to the last edge
    else
    {
        Edge *tmp = n1->list_edges;
        while (tmp->next != NULL)
            tmp = tmp->next;
        tmp->next = edge;
    }

    // Increment the number of edges of origin Node
    n1->num_edges++;
    // Increment the total number of edges of the graph
    total_edges++;
}

double calc_norm (const double x1, const double y1, const double z1,\
                  const double x2, const double y2, const double z2)
{
    return sqrt(pow(x2-x1,2) + pow(y2-y1,2) + pow(z2-z1,2));
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

void Graph::print ()
{
    	Node *ptr = list_nodes;
	printf("======================= PRINTING GRAPH ================================\n");
	while (ptr != NULL)
	{
        	printf("|| %d (%.2lf %.2lf %.2lf) d_ori(%.2lf %.2lf %.2lf) Edges = %d (%d) ||",ptr->id,ptr->x,ptr->y,ptr->z,\
											ptr->d_ori[0],ptr->d_ori[1],ptr->d_ori[2],\
											ptr->num_edges,ptr->is_terminal);
	
        	Edge *ptrl = ptr->list_edges;
		while (ptrl != NULL)
		{
			printf(" --> || %d %.2lf (%.2lf %.2lf %.2lf) ||",ptrl->id,ptrl->w,\
                                ptrl->dest->x,ptrl->dest->y,ptrl->dest->z);
			ptrl = ptrl->next;
		}
		printf("\n");
		ptr = ptr->next;
	}
	printf("=======================================================================\n");
    printf("Number of nodes = %d\n",total_nodes);
    printf("Number of edges = %d\n",total_edges);
}



