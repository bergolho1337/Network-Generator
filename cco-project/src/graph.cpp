#include "graph.h"

Graph::Graph ()
{
    last_node = NULL;
    list_nodes = NULL;
    total_nodes = 0;
    total_edges = 0;
}

/*
Graph::Graph (const char filename[])
{
    last_node = NULL;
    list_nodes = NULL;

    FILE *file = fopen(filename,"r");
    if (!file)
    {
        fprintf(stderr,"[-] ERROR! Cannot open filename '%s'\n",filename);
        exit(EXIT_FAILURE);
    }

    int num_nodes, num_edges;
    double pos[3];
    int dir[2];
    fscanf(file,"%d %d",&num_nodes,&num_edges);

    for (int i = 0; i < num_nodes; i++)
    {
        fscanf(file,"%lf %lf %lf",&pos[0],&pos[1],&pos[2]);
        insert_node_graph(pos);
    }
    for (int i = 0; i < num_edges; i++)
    {
        fscanf(file,"%d %d",&dir[0],&dir[1]);
        insert_edge_graph(dir[0],dir[1]);
    }

    // Unitary test
    assert(num_nodes == total_nodes);
    assert(num_edges == total_edges);

    fclose(file);
    
}
*/

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

void Graph::insert_edge_graph (const int id_1, const int id_2, const double Q_term, const double p_term)
{
    // TODO: Check if the edge is invalid
    //if (id_1 == id_2) return;

    Node *n1 = search_node(id_1);
    Node *n2 = search_node(id_2);

    double length = calc_norm(n1->x,n1->y,n1->z,n2->x,n2->y,n2->z);
    double radius = poiseuille(Q_term,p_term,length);
    //printf("Radius = %.10lf\n",radius);
    Edge *edge = new Edge(id_2,length,n2);

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
    fprintf(stderr,"[-] ERROR! Node %d not found!\n",id);

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

double calc_norm (const double x1, const double y1, const double z1,\
                  const double x2, const double y2, const double z2)
{
    return sqrt(pow(x2-x1,2) + pow(y2-y1,2) + pow(z2-z1,2));
}

// Return the radius of a segment with a prescribed Qterm, pterm and length
double poiseuille (const double Q_term, const double p_term, const double length)
{
    return pow( (8.0*Q_term*length*ETA)/(p_term*M_PI) , 0.25 );
}

/*

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