#include "cco.h"

struct cco_network* new_cco_network (struct user_options *options)
{
    struct cco_network *result = (struct cco_network*)malloc(sizeof(struct cco_network));
    result->point_list = new_point_list();
    result->segment_list = new_segment_list();
    result->num_terminals = 0;
    result->Q_perf = options->Q_perf;
    result->p_perf = options->p_perf;
    result->r_perf = options->r_perf;
    result->N_term = options->N_term;
    return result;
}

void grow_tree (struct cco_network *the_network)
{
    //test1(the_network);
    test2(the_network);
}

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
    double Q = the_network->Q_perf;
    double p = the_network->p_perf;

    struct point_list *p_list = the_network->point_list;
    struct segment_list *s_list = the_network->segment_list;

    // Root
    double pos1[3] = {0,0,0};
    double pos2[3] = {-3,-3,0};
    struct point_node *A = insert_point(p_list,pos1);
    struct point_node *B = insert_point(p_list,pos2);
    
    struct segment *iroot = new_segment(A,B,NULL,NULL,NULL,Q,p);
    struct segment_node *iroot_node = insert_segment_node(s_list,iroot);

    // First segment
    double pos3[3] = {1,-3,0};
    build_segment(the_network,0,pos3);

    // Second segment
    double pos4[3] = {-2,-4,0};
    build_segment(the_network,1,pos4);

    // Third segment
    double pos5[3] = {-1,-3,0};
    build_segment(the_network,2,pos5);

    print_list(p_list);
    print_list(s_list);
}

void build_segment (struct cco_network *the_network, const uint32_t index, const double new_pos[])
{
    struct point_list *p_list = the_network->point_list;
    struct segment_list *s_list = the_network->segment_list;
    double Q = the_network->Q_perf;
    double p = the_network->p_perf;

    struct segment_node *iconn = search_segment_node(s_list,index);

    // Create the middle point
    double middle_pos[3];
    calc_middle_point_segment(iconn,middle_pos);
    struct point_node *M = insert_point(p_list,middle_pos);

    // Create ibiff
    struct segment *ibiff = new_segment(M,iconn->value->dest,\
                            iconn->value->left,iconn->value->right,iconn,Q,p);
    struct segment_node *ibiff_node = insert_segment_node(s_list,ibiff);
    ibiff->ndist = iconn->value->ndist;

    // Create inew
    struct point_node *T = insert_point(p_list,new_pos);
    struct segment *inew = new_segment(M,T,\
                            NULL,NULL,iconn,Q,p);
    struct segment_node *inew_node = insert_segment_node(s_list,inew);
    
    // Update iconn pointers
    iconn->value->dest = M;
    iconn->value->left = ibiff_node;
    iconn->value->right = inew_node;

    // Update ndist
    struct segment_node *tmp = iconn;
    while (tmp != NULL)
    {
        tmp->value->ndist++;
        tmp = tmp->value->parent;
    }

}

void calc_middle_point_segment (struct segment_node *s, double pos[])
{
    struct point *src = s->value->src->value;
    struct point *dest = s->value->dest->value;

    pos[0] = (src->x + dest->x) / 2.0;
    pos[1] = (src->y + dest->y) / 2.0;
    pos[2] = (src->z + dest->z) / 2.0;
}

void write_to_vtk (struct cco_network *the_network)
{
    uint32_t num_points = the_network->point_list->num_nodes;
    struct point_list *p_list = the_network->point_list;
    struct point_node *p_tmp = p_list->list_nodes;
    uint32_t num_segments = the_network->segment_list->num_nodes;
    struct segment_list *s_list = the_network->segment_list;
    struct segment_node *s_tmp = s_list->list_nodes;

    FILE *file = fopen("cco_tree.vtk","w+");

    fprintf(file,"# vtk DataFile Version 3.0\n");
    fprintf(file,"Tree\n");
    fprintf(file,"ASCII\n");
    fprintf(file,"DATASET POLYDATA\n");
    fprintf(file,"POINTS %u float\n",num_points);
    while (p_tmp != NULL)
    {
        fprintf(file,"%g %g %g\n",p_tmp->value->x,p_tmp->value->y,p_tmp->value->z);
        p_tmp = p_tmp->next;
    }
        
    fprintf(file,"LINES %u %u\n",num_segments,num_segments*3);
    while (s_tmp != NULL)
    {
        fprintf(file,"2 %u %u\n",s_tmp->value->src->id,s_tmp->value->dest->id);
        s_tmp = s_tmp->next;
    }
    fclose(file);
}