#include "test.h"

void test1 (struct cco_network *the_network)
{
    double Q = the_network->Q_perf;
    double p = the_network->p_perf;

    struct point_list *p_list = the_network->point_list;

    // Positions
    double pos1[3] = {0,0,0};
    double pos2[3] = {-3,-3,0};
    double pos3[3] = {-1,-1,0};
    double pos4[3] = {1,-3,0};
    double pos5[3] = {-2,-2,0};
    double pos6[3] = {-2,-4,0};
    double pos7[3] = {0,-2,0};
    double pos8[3] = {-1,-3,0};

    struct point_node *A = insert_point(p_list,pos1);
    struct point_node *B = insert_point(p_list,pos2);
    struct point_node *C = insert_point(p_list,pos3);
    struct point_node *D = insert_point(p_list,pos4);
    struct point_node *E = insert_point(p_list,pos5);
    struct point_node *F = insert_point(p_list,pos6);
    struct point_node *G = insert_point(p_list,pos7);
    struct point_node *H = insert_point(p_list,pos8);

    print_list(p_list);

    // Segments
    struct segment_list *s_list = the_network->segment_list;

    struct segment *s1 = new_segment(A,C,NULL,NULL,NULL,Q,p);
    struct segment *s2 = new_segment(C,E,NULL,NULL,NULL,Q,p);
    struct segment *s3 = new_segment(C,G,NULL,NULL,NULL,Q,p);
    struct segment *s4 = new_segment(G,H,NULL,NULL,NULL,Q,p);
    struct segment *s5 = new_segment(G,D,NULL,NULL,NULL,Q,p);
    struct segment *s6 = new_segment(E,B,NULL,NULL,NULL,Q,p);
    struct segment *s7 = new_segment(E,F,NULL,NULL,NULL,Q,p);

    struct segment_node *s_node1 = insert_segment_node(s_list,s1);
    struct segment_node *s_node2 = insert_segment_node(s_list,s2);
    struct segment_node *s_node3 = insert_segment_node(s_list,s3);
    struct segment_node *s_node4 = insert_segment_node(s_list,s4);
    struct segment_node *s_node5 = insert_segment_node(s_list,s5);
    struct segment_node *s_node6 = insert_segment_node(s_list,s6);
    struct segment_node *s_node7 = insert_segment_node(s_list,s7);

    // Set pointers
    s1->parent = NULL;
    s1->left = s_node2;
    s1->right = s_node3;

    s2->parent = s_node1;
    s2->left = s_node6;
    s2->right = s_node7;

    s3->parent = s_node1;
    s3->left = s_node4;
    s3->right = s_node5;

    s4->parent = s_node3;
    s4->left = NULL;
    s4->right = NULL;

    s5->parent = s_node3;
    s5->left = NULL;
    s5->right = NULL;
    
    s6->parent = s_node2;
    s6->left = NULL;
    s6->right = NULL;

    s7->parent = s_node2;
    s7->left = NULL;
    s7->right = NULL;

    print_list(s_list);
}

void test2 (struct cco_network *the_network)
{
    double Q_perf = the_network->Q_perf;
    double p_perf = the_network->p_perf;
    double p_term = the_network->p_term;
    double delta_p = p_perf - p_term;

    struct point_list *p_list = the_network->point_list;
    struct segment_list *s_list = the_network->segment_list;

    // Root
    double pos1[3] = {0,0,0};
    double pos2[3] = {-3,-3,0};
    struct point_node *A = insert_point(p_list,pos1);
    struct point_node *B = insert_point(p_list,pos2);
    
    struct segment *iroot = new_segment(A,B,NULL,NULL,NULL,Q_perf,p_perf);
    struct segment_node *iroot_node = insert_segment_node(s_list,iroot);
    rescale_root(iroot_node,Q_perf,delta_p);
    the_network->num_terminals = 1;

    // First segment
    double pos3[3] = {1,-3,0};
    build_segment(the_network,NULL,0,pos3);

    // Second segment
    double pos4[3] = {-2,-4,0};
    build_segment(the_network,NULL,1,pos4);

    // Third segment
    double pos5[3] = {-1,-3,0};
    build_segment(the_network,NULL,2,pos5);

    print_list(p_list);
    print_list(s_list);
}

// Build the example network from Rafael's thesis
void test3 (struct cco_network *the_network)
{
    double Q_perf = the_network->Q_perf;
    double p_perf = the_network->p_perf;
    double p_term = the_network->p_term;
    double delta_p = p_perf - p_term;

    struct point_list *p_list = the_network->point_list;
    struct segment_list *s_list = the_network->segment_list;

    // Root
    double pos1[3] = {0,0,0};
    double pos2[3] = {0,-4,0};
    struct point_node *A = insert_point(p_list,pos1);
    struct point_node *B = insert_point(p_list,pos2);
    
    struct segment *iroot = new_segment(A,B,NULL,NULL,NULL,Q_perf,p_perf);
    struct segment_node *iroot_node = insert_segment_node(s_list,iroot);
    rescale_root(iroot_node,Q_perf,delta_p);
    the_network->num_terminals = 1;

    // First segment
    double pos3[3] = {2,-2.5,0};
    build_segment(the_network,NULL,0,pos3);

    // Second segment
    double pos4[3] = {-2,-1.5,0};
    build_segment(the_network,NULL,0,pos4);

    // Third segment
    double pos5[3] = {1,-5,0};
    build_segment(the_network,NULL,1,pos5);

    print_list(p_list);
    print_list(s_list);
}

void check_bifurcation_rule (struct cco_network *the_network)
{
    FILE *log_file = the_network->log_file;

    //printf("[!] Checking bifurcation rule!\n");
    fprintf(log_file,"[!] Checking bifurcation rule!\n");

    struct segment_list *s_list = the_network->segment_list;
    struct segment_node *tmp, *tmp_left, *tmp_right;
    
    tmp = s_list->list_nodes;
    while (tmp != NULL)
    {
        tmp_left = tmp->value->left;
        tmp_right = tmp->value->right;

        if (tmp_left && tmp_right)
        {
            double r = pow(tmp->value->radius,GAMMA);
            double r_left = pow(tmp_left->value->radius,GAMMA);
            double r_right = pow(tmp_right->value->radius,GAMMA);

            double diff = r_left + r_right - r;
            if (diff > TOLERANCE)
            {
                fprintf(stderr,"[test] Error when checking bifurcation rule !\n");
            }
            fprintf(log_file,"%g = %g + %g --> %g = %g\n",r,r_left,r_right,r,r_left + r_right);
        }
        tmp = tmp->next;
    }
}