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
    test1(the_network);
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

    struct segment_list *s_list = the_network->segment_list;

    struct segment *s1 = new_segment(A,C,NULL,NULL,NULL,Q,p);
    struct segment *s2 = new_segment(C,E,NULL,NULL,NULL,Q,p);
    struct segment *s3 = new_segment(C,G,NULL,NULL,NULL,Q,p);
    struct segment *s4 = new_segment(G,H,NULL,NULL,NULL,Q,p);
    struct segment *s5 = new_segment(G,D,NULL,NULL,NULL,Q,p);
    struct segment *s6 = new_segment(E,B,NULL,NULL,NULL,Q,p);
    struct segment *s7 = new_segment(E,F,NULL,NULL,NULL,Q,p);

    insert_segment_node(s_list,s1);
    insert_segment_node(s_list,s2);
    insert_segment_node(s_list,s3);
    insert_segment_node(s_list,s4);
    insert_segment_node(s_list,s5);
    insert_segment_node(s_list,s6);
    insert_segment_node(s_list,s7);

    print_list(s_list);

    /*
    Segment s1(&A,&C,NIL,NIL,NIL,Q_perf,p_perf);
    Segment s2(&C,&E,NIL,NIL,NIL,Q_perf,p_perf);
    Segment s3(&C,&G,NIL,NIL,NIL,Q_perf,p_perf);
    Segment s4(&G,&H,NIL,NIL,NIL,Q_perf,p_perf);
    Segment s5(&G,&D,NIL,NIL,NIL,Q_perf,p_perf);
    Segment s6(&E,&B,NIL,NIL,NIL,Q_perf,p_perf);
    Segment s7(&E,&F,NIL,NIL,NIL,Q_perf,p_perf);

    segments.push_back(s1);
    segments.push_back(s2);
    segments.push_back(s3);
    segments.push_back(s4);
    segments.push_back(s5);
    segments.push_back(s6);
    segments.push_back(s7);
    */
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