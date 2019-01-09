#include "../include/cco.h"

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

void make_root (Graph *the_network, const double Q_perf, const double Nterm)
{
    double rsupp = sqrt(Q_perf / (Nterm * M_PI));
    double xdist, ydist, rand_number;
    
    srand(time(NULL));
    do
    {
        rand_number = (double)rand() / (double)RAND_MAX;
        xdist = -rsupp + rand_number*2.0*rsupp;

        rand_number = (double)rand() / (double)RAND_MAX;
        ydist = 0.0 + rand_number*(-2.0)*rsupp;
         
    }while (!is_inside_circle(xdist,ydist,rsupp));

    insert_node_graph(the_network,0.0,0.0);
    insert_node_graph(the_network,xdist,ydist);

    insert_edge_graph(the_network,0,1);

}

void grow_cco_tree (Graph *the_network, User_Options *options)
{
    make_root(the_network,options->Q_perf,options->N_term);
}