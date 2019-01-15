#include "../include/cco.h"

CCO_Network::CCO_Network (User_Options *options)
{
    num_terminals = 0;
    Q_perf = options->Q_perf;
    p_perf = options->p_perf;
    r_perf = options->r_perf;
}

void CCO_Network::grow_tree ()
{
    make_root();
}

void CCO_Network::make_root ()
{
    
    // Calculating the radius of the first microcirculatory black-box (Nterm = 1 -> root)
    double r_supp = sqrt(Q_perf / M_PI);

    Point proximal; proximal.x = 0.0; proximal.y = 0.0; proximal.z = 0.0; 
    Point distal;
    
    // Calculate teh distal position of the root inside the circle with radius r_supp
    srand(time(NULL));
    generate_point_inside_circle(&distal,r_supp);
    
    // Insert the root into the graph
    //the_network->insert_node_graph(pos_prox);
    //the_network->insert_node_graph(pos_dist);

    // (Nterm = 1) in this case
    //the_network->insert_edge_graph(0,1,Q_perf,p_perf);
}

bool is_inside_circle (Point *p, const Point c, const double radius)
{
    double d = sqrt(pow(p->x - c.x,2) + pow(p->y - c.y,2) + pow(p->z - c.z,2));

    return (d <= radius) ? true : false;
}

void generate_point_inside_circle (Point *p, const double radius)
{
    // Center of the perfusion circle
    Point center; 
    center.x = 0.0;
    center.y = -radius;
    center.z = 0.0;

    double rand_number;
    do
    {
        rand_number = (double)rand() / (double)RAND_MAX;
        p->x = -radius + 2.0*rand_number*radius;

        rand_number = (double)rand() / (double)RAND_MAX;
        p->y = -2.0*rand_number*radius;

        p->z = 0.0;
         
    }while (!is_inside_circle(p,center,radius));
}

/*
bool is_inside_circle (const double x, const double y, const double r)
{
    // Center of the perfusion circle
    double center_x = 0.0;
    double center_y = -r;

    double d = sqrt(pow(x-center_x,2) + pow(y-center_y,2));

    printf("x = %lf\n",x);
    printf("y = %lf\n",y);
    printf("center_x = %lf\n",center_x);
    printf("center_y = %lf\n\n",center_y);

    return (d <= r) ? true : false;
}

// This function implements the part:
// "A. Planting the Root" from the article
void make_root (Graph *the_network, const double Q_perf, const double p_perf)
{
    // Calculating the radius of the first microcirculatory black-box (Nterm = 1 -> root)
    double rsupp = sqrt(Q_perf / (1 * M_PI));

    double rand_number;
    double pos_prox[2] = {0,0}; 
    double pos_dist[2];
    
    // Calculate teh distal position of the root inside the circle with radius r_supp
    srand(time(NULL));
    do
    {
        rand_number = (double)rand() / (double)RAND_MAX;
        pos_dist[0] = -rsupp + rand_number*2.0*rsupp;

        rand_number = (double)rand() / (double)RAND_MAX;
        pos_dist[1] = 0.0 + rand_number*(-2.0)*rsupp;
         
    }while (!is_inside_circle(pos_dist[0],pos_dist[1],rsupp));

    // Insert the root into the graph
    the_network->insert_node_graph(pos_prox);
    the_network->insert_node_graph(pos_dist);

    // (Nterm = 1) in this case
    the_network->insert_edge_graph(0,1,Q_perf,p_perf);

}

// Construir um vetor de Segments
void add_terminal (Graph *the_network, const int num_term)
{
    inflate_support_circle(num_term);

    // Stretch the coordinates of the others points
    streatch_coordinate_points(the_network);

    // Choose the position of the new terminal (must be inside the new supporting circle)
    calculate_new_terminal_position();

    double d_proj, d_crit, eval;
    double min_eval;
    int min_segment_index;
    
    // For each segment, test if it will be the best one
    for (int i = 0; i < num_segments; i++)
    {
        d_proj = segment[i].calc_dproj(x,y);

        if (d_proj >= 0 && d_proj <= 1)
            d_crit = segment[i].calc_dortho(x,y);
        else
            d_crit = segment[i].calc_dend(x,y);

        if (d_crit > d_threash)
            eval = build_segment(x,y,segment);
            if (eval < min_eval)
                min_eval = eval;
                min_segment_index = i;
            destroy_segment(segment);    

    }
    eval = build_segment(x,y,segment[min_segment_index]);
}

void generate_terminals (Graph *the_network, const int N_term)
{
    int num_term = 1;
    while (num_term <= N_term)
    {
        printf("[+] Growing terminal %d\n",num_term);
        
        add_terminal(the_network,num_term);
        
        num_term++;
    }
}

void grow_cco_tree (Graph *the_network, User_Options *options)
{
    make_root(the_network,options->Q_perf,options->p_perf);

    generate_terminals(the_network,options->N_term);
}
*/