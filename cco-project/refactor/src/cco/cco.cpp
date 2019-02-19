#include "cco.h"

struct cco_network* new_cco_network (struct user_options *options)
{
    struct cco_network *result = (struct cco_network*)malloc(sizeof(struct cco_network));
    result->point_list = new_point_list();
    result->segment_list = new_segment_list();
    result->num_terminals = 0;
    result->Q_perf = options->Q_perf;
    result->p_perf = options->p_perf;
    result->p_term = options->p_term;
    result->r_perf = options->r_perf;
    result->N_term = options->N_term;
    return result;
}

void grow_tree (struct cco_network *the_network)
{
    //test1(the_network);
    //test2(the_network);
    //test3(the_network);
    test_cco(the_network);

    check_bifurcation_rule(the_network);
}

void check_bifurcation_rule (struct cco_network *the_network)
{
    printf("[!] Checking bifurcation rule!\n");

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

            printf("%g = %g + %g --> %g = %g\n",r,r_left,r_right,r,r_left + r_right);
        }
        tmp = tmp->next;
    }
}

bool collision_detection (const double x1, const double y1, const double z1,\
                          const double x2, const double y2, const double z2,\
                          const double x3, const double y3, const double z3,\
                          const double x4, const double y4, const double z4)
{
    double denominator = ((x2 - x1) * (y3 - y4)) - ((y2 - y1) * (x3 - x4));
    double numerator1 = ((y1 - y4) * (x3 - x4)) - ((x1 - x4) * (y3 - y4));
    double numerator2 = ((y1 - y4) * (x2 - x1)) - ((x1 - x4) * (y2 - y1));
    
    // Detect coincident lines (has a problem, read below)
    if (denominator == 0) return numerator1 == 0 && numerator2 == 0;

    double r = numerator1 / denominator;
    double s = numerator2 / denominator;

    return (r >= 0 && r <= 1) && (s >= 0 && s <= 1);
}

bool has_collision (struct segment_list *s_list, struct segment_node *s, const double new_pos[])
{
    double middle_pos[3];
    calc_middle_point_segment(s,middle_pos); 

    printf("[+] Trying connection with segment %u\n",s->id);
    struct segment_node *tmp = s_list->list_nodes;
    while (tmp != NULL)
    {
        if (tmp->id != s->id)
        {
            printf("\t[!] Checking collison between segment %d\n",tmp->id);

            // Get the reference to the points of the current segment
            struct point *src = tmp->value->src->value;
            struct point *dest = tmp->value->dest->value;

            bool intersect = collision_detection(middle_pos[0],middle_pos[1],middle_pos[2],\
                                            new_pos[0],new_pos[1],new_pos[2],\
                                            src->x,src->y,src->z,\
                                            dest->x,dest->y,dest->z);  

            if (intersect)
            {
                printf("\t[-] ERROR! Intersection with segment %d !\n",tmp->id);
                return true;
            }          
        }
        tmp = tmp->next;
    }
    return false;
}

bool check_collisions (struct cco_network *the_network, const double new_pos[])
{
    printf("[!] Checking collisions!\n");

    struct segment_list *s_list = the_network->segment_list;
    struct segment_node *tmp = s_list->list_nodes;

    while (tmp != NULL)
    {
        if (has_collision(s_list,tmp,new_pos))
            return false;

        tmp = tmp->next;
    }
    return true;
}

void rescale_root (struct segment_node *iroot, const double Q_perf, const double delta_p)
{
    // Set the flux and pressure of this segment
    iroot->value->Q = Q_perf;
    iroot->value->p = delta_p;

    // Calculate the segment size
    struct point *prox = iroot->value->src->value;
    struct point *dist = iroot->value->dest->value;
    double length = euclidean_norm(prox->x,prox->y,prox->z,dist->x,dist->y,dist->z);
    iroot->value->length = length;

    // Calculate the segment relative resistance
    double R = 8.0 * ETA * length / M_PI;
    iroot->value->resistance = R;

    // Calculate the radius using Poisselle's law
    iroot->value->radius = pow(R * Q_perf / delta_p , 0.25);

}

void rescale_tree (struct segment_node *ibiff, struct segment_node *iconn, struct segment_node *inew,\
                 const double Q_perf, const double delta_p, const int num_terminals)
{

    double Q_term = Q_perf / num_terminals;
    
    // inew: Calculate resistance, flux and radius using (2.4) from Rafael's thesis
    calc_relative_resistance_term(inew);
    calc_radius_term(inew,Q_term,delta_p);
    
    // iconn + inew: Calculate bifurcation ratio using (2.31)
    inew->value->beta = calc_bifurcation_ratio(iconn->value->radius,inew->value->radius,true);
    iconn->value->beta = calc_bifurcation_ratio(iconn->value->radius,inew->value->radius,false);

    // TODO: Need to revise this ...
    // Preciso recalcular tudo do iconn ????
    // iconn: Recalculate resistance with the new bifurcation ratio
    if (iconn->value->left == NULL && iconn->value->right == NULL)
    {
        calc_relative_resistance_term(iconn);
        calc_radius_term(iconn,Q_term,delta_p);
    }

    // ibiff: Calculate resistance using (2.5), flux and radius using (2.4) from Rafael's thesis
    calc_relative_resistance_subtree(ibiff,iconn,inew);
    calc_radius_term(ibiff,Q_term,delta_p);

    // Rescale the until we reach the root by using the "parent" pointer
    struct segment_node *ipar = ibiff->value->parent;
    if (ipar != NULL)
    {
        struct segment_node *ipar_left = ipar->value->left;
        struct segment_node *ipar_right = ipar->value->right;
        rescale_until_root(ipar,ipar_left,ipar_right,Q_perf,delta_p,num_terminals);
    }
    

}

void rescale_until_root (struct segment_node *ipar, struct segment_node *ipar_left, struct segment_node *ipar_right,\
                        const double Q_perf, const double delta_p, const int num_terminals)
{
    printf("On rescale_until_root()\n");

    // Reach the root
    if (ipar == NULL) return;

    double Q_term = Q_perf / num_terminals;
    if (ipar_left != NULL && ipar_right != NULL)
    {
        
        // Recalculate bifurcation ratios for the offsprings using (2.31)
        ipar_right->value->beta = calc_bifurcation_ratio(ipar_left->value->radius,ipar_right->value->radius,true);
        ipar_left->value->beta = calc_bifurcation_ratio(ipar_left->value->radius,ipar_right->value->radius,false);

        // Recalculate resistance using (2.5)
        calc_relative_resistance_subtree(ipar,ipar_left,ipar_right);

        // Recalculate radius using (2.4)
        calc_radius_term(ipar,Q_term,delta_p);

        // Call the function recursively until we reach the root
        if (ipar->value->parent != NULL) 
            rescale_until_root(ipar->value->parent,\
                               ipar->value->parent->value->left,\
                               ipar->value->parent->value->right,\
                               Q_perf,delta_p,num_terminals); 
    }

}

void build_segment (struct cco_network *the_network, const uint32_t index, const double new_pos[])
{
    
    struct point_list *p_list = the_network->point_list;
    struct segment_list *s_list = the_network->segment_list;
    double Q_perf = the_network->Q_perf;
    double p_perf = the_network->p_perf;
    double p_term = the_network->p_term;
    double delta_p = p_perf - p_term;

    struct segment_node *iconn_node = search_segment_node(s_list,index);
    struct segment *iconn = iconn_node->value;

    // Create the middle point
    double middle_pos[3];
    calc_middle_point_segment(iconn_node,middle_pos);
    struct point_node *M = insert_point(p_list,middle_pos);

    // Create ibiff
    struct segment *ibiff = new_segment(M,iconn_node->value->dest,\
                            iconn_node->value->left,iconn_node->value->right,iconn_node,Q_perf,p_perf);
    struct segment_node *ibiff_node = insert_segment_node(s_list,ibiff);
    ibiff->ndist = iconn_node->value->ndist;
    ibiff->radius = iconn->radius;          // Copy from iconn !!!
    ibiff->resistance = iconn->resistance;  // Copy from iconn !!!

    // Create inew
    struct point_node *T = insert_point(p_list,new_pos);
    struct segment *inew = new_segment(M,T,\
                            NULL,NULL,iconn_node,Q_perf,p_perf);
    struct segment_node *inew_node = insert_segment_node(s_list,inew);
    
    // Update iconn pointers
    iconn_node->value->dest = M;
    iconn_node->value->left = ibiff_node;   // CONVENTION: Left will point to subtree
    iconn_node->value->right = inew_node;   // CONVENTION: Right will point to terminal

    // Update ndist
    struct segment_node *tmp = iconn_node;
    while (tmp->value->parent != NULL)
    {
        tmp->value->ndist++;
        tmp = tmp->value->parent;
    }
    tmp->value->ndist++;
    the_network->num_terminals = tmp->value->ndist;

    // "iconn" is equivalent to "ibiff" on Rafael's thesis and
    // "ibiff" is equivalent to "iconn"
    rescale_tree(iconn_node,ibiff_node,inew_node,Q_perf,delta_p,the_network->num_terminals);
}

void calc_relative_resistance_subtree (struct segment_node *ibiff, struct segment_node *iconn, struct segment_node *inew)
{
    calc_relative_resistance_term(ibiff);
    double R = ibiff->value->resistance;
    double R_left = pow(iconn->value->beta,4) / iconn->value->resistance;
    double R_right = pow(inew->value->beta,4) / inew->value->resistance;

    ibiff->value->resistance = R + pow( R_left + R_right , -1.0 );
}

void calc_relative_resistance_term (struct segment_node *iterm)
{
    struct point *src = iterm->value->src->value;
    struct point *dest = iterm->value->dest->value;
    double length = euclidean_norm(src->x,src->y,src->z,dest->x,dest->y,dest->z);

    iterm->value->resistance = 8.0 * ETA * length / M_PI;
}

void calc_radius_term (struct segment_node *iterm, const double Q_term, const double delta_p)
{
    double R = iterm->value->resistance;
    double Q = iterm->value->ndist * Q_term;

    iterm->value->radius = pow( R * Q / delta_p , 0.25 );
}

void calc_middle_point_segment (struct segment_node *s, double pos[])
{
    struct point *src = s->value->src->value;
    struct point *dest = s->value->dest->value;

    pos[0] = (src->x + dest->x) / 2.0;
    pos[1] = (src->y + dest->y) / 2.0;
    pos[2] = (src->z + dest->z) / 2.0;
}

double calc_bifurcation_ratio (const double r_left, const double r_right, const bool sign)
{
    // TODO: Calcular uma unica vez
    double expoent = -1.0/GAMMA;
    double base = r_right / r_left;

    if (sign)
        return pow( 1.0 + ( pow( base, -GAMMA) ) , expoent);
    else
        return pow( 1.0 + ( pow( base, GAMMA) ) , expoent);
}

void write_to_vtk (struct cco_network *the_network)
{
    uint32_t num_points = the_network->point_list->num_nodes;
    struct point_list *p_list = the_network->point_list;
    struct point_node *p_tmp = p_list->list_nodes;
    uint32_t num_segments = the_network->segment_list->num_nodes;
    struct segment_list *s_list = the_network->segment_list;
    struct segment_node *s_tmp = s_list->list_nodes;

    FILE *file = fopen("output/cco_tree.vtk","w+");

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
    fprintf(file,"CELL_DATA %u\n",num_segments);
    fprintf(file,"SCALARS radius float\n");
    fprintf(file,"LOOKUP_TABLE default\n");
    s_tmp = s_list->list_nodes;
    while (s_tmp != NULL)
    {
        fprintf(file,"%g\n",s_tmp->value->radius);
        s_tmp = s_tmp->next;
    }
    fclose(file);
}

void draw_perfusion_area (const double radius)
{
    // Create a circle
    vtkSmartPointer<vtkRegularPolygonSource> polygonSource =
      vtkSmartPointer<vtkRegularPolygonSource>::New();
    
    polygonSource->GeneratePolygonOff(); // Outline of the circle
    polygonSource->SetNumberOfSides(50);
    polygonSource->SetRadius(radius);
    polygonSource->SetCenter(0,-radius,0);

    vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
    writer->SetFileName("output/perfusion_area.vtp");
    writer->SetInputConnection(polygonSource->GetOutputPort());
    writer->Write();
}

// --------------------------------------------------------------------------------------------
// TEST FUNCTIONS
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
    build_segment(the_network,0,pos3);

    // Second segment
    double pos4[3] = {-2,-1.5,0};
    build_segment(the_network,0,pos4);

    // Third segment
    double pos5[3] = {1,-5,0};
    build_segment(the_network,1,pos5);

    print_list(p_list);
    print_list(s_list);
}

void make_root (struct cco_network *the_network, const double r_supp)
{
    double Q_perf = the_network->Q_perf;
    double p_perf = the_network->p_perf;
    double p_term = the_network->p_term;
    double r_perf = the_network->r_perf;
    double delta_p = p_perf - p_term;

    struct point_list *p_list = the_network->point_list;
    struct segment_list *s_list = the_network->segment_list;

    // Positions from the root
    double x_inew[3] = {0,0,0};
    double x_prox[3] = {0,0,0};

    // Sort the distal position of the root until its size is larger than the perfusion radius  
    while (euclidean_norm(x_prox[0],x_prox[1],x_prox[2],x_inew[0],x_inew[1],x_inew[2]) < r_perf)
        generate_point_inside_perfusion_area(x_inew,r_supp);

    // Insert points and create the root segment
    struct point_node *A = insert_point(p_list,x_prox);
    struct point_node *B = insert_point(p_list,x_inew);

    struct segment *iroot = new_segment(A,B,NULL,NULL,NULL,Q_perf,p_perf);
    struct segment_node *iroot_node = insert_segment_node(s_list,iroot);
    rescale_root(iroot_node,Q_perf,delta_p);
    the_network->num_terminals = 1;
}

void generate_terminal (struct cco_network *the_network)
{
    int num_terminals = the_network->num_terminals;
    double r_perf = the_network->r_perf;

    double new_pos[3];
    bool point_is_ok = false;
    uint32_t tosses = 0;
    double d_threash = calc_dthreashold(r_perf,num_terminals);
    
    while (!point_is_ok)
    {
        // Generate the terminal position inside the perfusion area
        generate_point_inside_perfusion_area(new_pos,r_perf);

        // Check the distance criterion for this point
        //ret = connection_search_paper(pos,d_threash);

        // Check collision with other segments
        point_is_ok = check_collisions(the_network,new_pos);

        // If the point does not attend the distance criterion or if there is a collision
        // we need to choose another point.
        if (!point_is_ok)
        {
            tosses++;
            if (tosses > NTOSS)
            {
                printf("[!] Reducing dthreash! Before = %.2lf || Now = %.2lf \n",\
                        d_threash,d_threash*0.9);
                d_threash *= 0.9;
                tosses = 0;
            }
        }
    }
}

void test_cco (struct cco_network *the_network)
{
    double Q_perf = the_network->Q_perf;
    double p_perf = the_network->p_perf;
    double p_term = the_network->p_term;
    double r_perf = the_network->r_perf;
    double delta_p = p_perf - p_term;        

    struct point_list *p_list = the_network->point_list;
    struct segment_list *s_list = the_network->segment_list;

    draw_perfusion_area(r_perf);

    make_root(the_network,r_perf);

    // Main iteration loop
    while (the_network->num_terminals <= the_network->N_term)
    {
        printf("%s\n",PRINT_LINE);
        printf("[!] Working on terminal number %d\n",the_network->num_terminals);            

        generate_terminal(the_network);

        the_network->num_terminals++;

        printf("%s\n",PRINT_LINE);
    }

    print_list(p_list);
    print_list(s_list);


}

double calc_dthreashold (const double radius, const int num_terminals)
{
    return sqrt( (double)(M_PI*radius*radius) / (double)(num_terminals) );
}