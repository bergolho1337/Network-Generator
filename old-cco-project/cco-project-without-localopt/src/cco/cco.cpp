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
    result->A_perf = M_PI * result->r_perf * result->r_perf;
    result->log_file = fopen("output.log","w+");

    uint32_t size = strlen(options->config->name) + 1;
    result->cost_function_name = (char*)malloc(sizeof(char)*size);
    strcpy(result->cost_function_name,options->config->name);
    printf("[cco] Cost function name = %s\n",result->cost_function_name);

    if (options->use_cloud_points)
    {
        result->using_cloud_points = true;
        size = strlen(options->cloud_points_filename) + 1;
        result->cloud_points_filename = (char*)malloc(sizeof(char)*size);
        strcpy(result->cloud_points_filename,options->cloud_points_filename);

        printf("[cco] Using cloud of points\n");
        printf("[cco] Cloud points filename :> \"%s\"\n",result->cloud_points_filename);
    }
    else
    {
        result->using_cloud_points = false;
        
        printf("[cco] No cloud of points provided\n");
        printf("[cco] Generating cloud of points based on value of \"r_perf\"\n");    
    }

    return result;
}

void free_cco_network (struct cco_network *the_network)
{
    free_point_list(the_network->point_list);
    free_segment_list(the_network->segment_list);
    free(the_network->cost_function_name);
    if (the_network->using_cloud_points)
        free(the_network->cloud_points_filename);
    fclose(the_network->log_file);

    free(the_network);
}

void grow_tree (struct cco_network *the_network, struct user_options *options)
{
    printf("\n[cco] Growing CCO network !\n");

    if (!the_network->using_cloud_points)
        grow_tree_default(the_network,options);
    else
        grow_tree_using_cloud_points(the_network,options);
    
    // Unitary test
    check_bifurcation_rule(the_network);

}

void grow_tree_default (struct cco_network *the_network, struct user_options *options)
{
    FILE *log_file = the_network->log_file;

    struct cost_function_config *config = options->config;

    double Q_perf = the_network->Q_perf;
    double p_perf = the_network->p_perf;
    double p_term = the_network->p_term;
    double r_perf = the_network->r_perf;
    double delta_p = p_perf - p_term;        

    struct point_list *p_list = the_network->point_list;
    struct segment_list *s_list = the_network->segment_list;

    make_root(the_network);

    // Main iteration loop
    while (the_network->num_terminals < the_network->N_term)
    {
        printf("%s\n",PRINT_LINE);
        printf("[cco] Working on terminal number %d\n",the_network->num_terminals);            
        fprintf(log_file,"%s\n",PRINT_LINE);
        fprintf(log_file,"[cco] Working on terminal number %d\n",the_network->num_terminals);

        generate_terminal(the_network,config);

        printf("%s\n",PRINT_LINE);
        fprintf(log_file,"%s\n",PRINT_LINE);
    }

    draw_perfusion_area(the_network);

    //print_list(p_list);
    //print_list(s_list);

    // Write to the logfile
    fprintf(log_file,"%s\n",PRINT_LINE);
    write_list(p_list,log_file);
    fprintf(log_file,"%s\n",PRINT_LINE);

    fprintf(log_file,"%s\n",PRINT_LINE);
    write_list(s_list,log_file);
    fprintf(log_file,"%s\n",PRINT_LINE);
}

void grow_tree_using_cloud_points (struct cco_network *the_network, struct user_options *options)
{
    FILE *log_file = the_network->log_file;

    struct cost_function_config *config = options->config;

    double Q_perf = the_network->Q_perf;
    double p_perf = the_network->p_perf;
    double p_term = the_network->p_term;
    double r_perf = the_network->r_perf;
    double delta_p = p_perf - p_term;        

    struct point_list *p_list = the_network->point_list;
    struct segment_list *s_list = the_network->segment_list;

    if (the_network->using_cloud_points)
    {
        std::vector<struct point> cloud_points; 
        read_cloud_points(the_network->cloud_points_filename,cloud_points);

        make_root_using_cloud_points(the_network,cloud_points);

        // Main iteration loop
        while (the_network->num_terminals < the_network->N_term)
        {
            printf("%s\n",PRINT_LINE);
            printf("[cco] Working on terminal number %d\n",the_network->num_terminals);            
            fprintf(log_file,"%s\n",PRINT_LINE);
            fprintf(log_file,"[cco] Working on terminal number %d\n",the_network->num_terminals);

            generate_terminal_using_cloud_points(the_network,config,cloud_points);

            printf("%s\n",PRINT_LINE);
            fprintf(log_file,"%s\n",PRINT_LINE);
        }
    }

    //print_list(p_list);
    //print_list(s_list);

    // Write to the logfile
    fprintf(log_file,"%s\n",PRINT_LINE);
    write_list(p_list,log_file);
    fprintf(log_file,"%s\n",PRINT_LINE);

    fprintf(log_file,"%s\n",PRINT_LINE);
    write_list(s_list,log_file);
    fprintf(log_file,"%s\n",PRINT_LINE);
}

bool has_collision (struct segment_list *s_list, struct segment_node *s, const double new_pos[], FILE *log_file)
{
    double middle_pos[3];
    calc_middle_point_segment(s,middle_pos); 

    //printf("[+] Trying connection with segment %u\n",s->id);
    //fprintf(log_file,"[+] Trying connection with segment %u\n",s->id);

    struct segment_node *tmp = s_list->list_nodes;
    while (tmp != NULL)
    {
        if (tmp->id != s->id)
        {
            //printf("\t[!] Checking collison between segment %d\n",tmp->id);
            //fprintf(log_file,"\t[!] Checking collison between segment %d\n",tmp->id);

            // Get the reference to the points of the current segment
            struct point *src = tmp->value->src->value;
            struct point *dest = tmp->value->dest->value;

            bool intersect = collision_detection(middle_pos[0],middle_pos[1],middle_pos[2],\
                                            new_pos[0],new_pos[1],new_pos[2],\
                                            src->x,src->y,src->z,\
                                            dest->x,dest->y,dest->z);  

            if (intersect)
            {
                //printf("\t[-] ERROR! Intersection with segment %d !\n",tmp->id);
                //fprintf(log_file,"\t[-] ERROR! Intersection with segment %d !\n",tmp->id);
                return true;
            }          
        }
        tmp = tmp->next;
    }
    return false;
}

bool check_collisions (struct cco_network *the_network, const double new_pos[],\
                    std::vector<struct segment_node*> &feasible_segments)
{
    FILE *log_file = the_network->log_file;

    //printf("[!] Checking collisions!\n");
    fprintf(log_file,"[!] Checking collisions!\n");

    struct segment_list *s_list = the_network->segment_list;
    struct segment_node *tmp = s_list->list_nodes;

    while (tmp != NULL)
    {
        if (!has_collision(s_list,tmp,new_pos,log_file))
        {
            feasible_segments.push_back(tmp);
        }

        tmp = tmp->next;
    }

    if (feasible_segments.size() == 0)
        return false;
    else
        return true;
}

bool connection_search (struct cco_network *the_network, const double pos[], const double d_threash)
{
    struct segment_list *s_list = the_network->segment_list;
    struct segment_node *tmp = s_list->list_nodes;

    while (tmp != NULL)
    {
        if (!distance_criterion(tmp,pos,d_threash))
            return false;

        tmp = tmp->next;
    }

    return true;
}

bool distance_criterion (struct segment_node *s, const double pos[], const double d_threash)
{
    double d_proj = calc_dproj(s,pos);
    
    double d_crit;
    if (d_proj >= 0 && d_proj <= 1)
        d_crit = calc_dortho(s,pos);
    else
        d_crit = calc_dend(s,pos);
    
    return (d_crit > d_threash) ? true : false;
}

bool has_deviation (struct segment_list *s_list, struct segment_node *inew,\
                    const double new_at, const double limit,\
                    const double c, const double cm, const double rc, const double rm)
{
    struct segment_node *tmp = s_list->list_nodes;

    while (tmp != NULL)
    {
        if (tmp != inew && is_terminal(tmp))
        {
            double at = calc_terminal_activation_time(tmp,c,cm,rc,rm);
            if (fabs(new_at - at) > limit)
                return true;
        }
        tmp = tmp->next;
    }

    return false;
}

bool is_terminal (struct segment_node *s)
{
    if (s->value->left == NULL && s->value->right == NULL)
        return true;
    else
        return false;
}

void rescale_root (struct segment_node *iroot, const double Q_perf, const double delta_p)
{
    // Set the flux and pressure of this segment
    iroot->value->Q = Q_perf;
    iroot->value->delta_p = delta_p;

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
    
    // inew: Calculate resistance, pressure, flux and radius 
    // using (2.32) and (2.22) from Rafael's thesis
    calc_relative_resistance_term(inew);
    calc_pressure_drop_term(inew,Q_term);

    // iconn: Recalculate resistance using (2.5) and pressure drop (2.7)
    if (iconn->value->left == NULL && iconn->value->right == NULL)
    {
        calc_relative_resistance_term(iconn);
        calc_pressure_drop_term(iconn,Q_term);
    }
    else
    {
        calc_relative_resistance_subtree(iconn,iconn->value->right,iconn->value->left);
        calc_pressure_drop_subtree(iconn,Q_term);
    }
        
    // Compute the radius ratio between left and right subtree from ibiff (2.32)
    double radius_ratio = calc_radius_ratio(iconn,inew,Q_term);

    // iconn + inew: Calculate bifurcation ratio using (2.31)
    inew->value->beta = calc_bifurcation_ratio(radius_ratio,false);
    iconn->value->beta = calc_bifurcation_ratio(radius_ratio,true);

    // ibiff: Calculate resistance using (2.5) and pressure drop using (2.7)
    calc_relative_resistance_subtree(ibiff,iconn,inew);
    calc_pressure_drop_subtree(ibiff,Q_term);

    // Rescale the until we reach the root by using the "parent" pointer
    struct segment_node *ipar = ibiff->value->parent;
    if (ipar != NULL)
    {
        struct segment_node *ipar_left = ipar->value->left;
        struct segment_node *ipar_right = ipar->value->right;
        rescale_until_root(ipar,ipar_left,ipar_right,Q_perf,delta_p,num_terminals);
    }
    // We are already at the root, so recalculate the radius using (2.19)
    else
    {
        ibiff->value->radius = pow(ibiff->value->resistance * Q_perf / delta_p , 0.25);
    }
}

void rescale_until_root (struct segment_node *ipar, struct segment_node *ipar_left, struct segment_node *ipar_right,\
                        const double Q_perf, const double delta_p, const int num_terminals)
{

    // Reach the root
    if (ipar == NULL) return;

    double Q_term = Q_perf / num_terminals;
    
    if (ipar_left != NULL && ipar_right != NULL)
    {

        // Calculate radius ratio between left and right subtree
        double radius_ratio = calc_radius_ratio(ipar_right,ipar_left,Q_term);

        // Recalculate bifurcation ratios for the offsprings using (2.31)
        ipar_left->value->beta = calc_bifurcation_ratio(radius_ratio,false);
        ipar_right->value->beta = calc_bifurcation_ratio(radius_ratio,true);

        // Recalculate resistance using (2.5) and pressure drop using (2.7)
        calc_relative_resistance_subtree(ipar,ipar_left,ipar_right);
        calc_pressure_drop_subtree(ipar,Q_term);

        // Call the function recursively until we reach the root
        if (ipar->value->parent != NULL) 
            rescale_until_root(ipar->value->parent,\
                               ipar->value->parent->value->left,\
                               ipar->value->parent->value->right,\
                               Q_perf,delta_p,num_terminals);
        // Recalculate the root radius when we reach this segment using (2.19)
        else
            ipar->value->radius = pow(ipar->value->resistance * Q_perf / delta_p , 0.25); 
    }

}

struct segment_node* build_segment (struct cco_network *the_network, const uint32_t index, const double new_pos[])
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
    struct segment *ibiff = new_segment(iconn_node->value->src,M,\
                            NULL,NULL,iconn_node->value->parent,Q_perf,p_perf);
    struct segment_node *ibiff_node = insert_segment_node(s_list,ibiff);
    ibiff->ndist = iconn_node->value->ndist;

    // Create inew
    struct point_node *T = insert_point(p_list,new_pos);
    struct segment *inew = new_segment(M,T,\
                            NULL,NULL,ibiff_node,Q_perf,p_perf);
    struct segment_node *inew_node = insert_segment_node(s_list,inew);
    
    // Update pointers
    //  iconn:
    iconn_node->value->parent = ibiff_node;
    iconn_node->value->src = M;

    //  ibiff:
    ibiff_node->value->left = inew_node;     // CONVENTION: Left will point to terminal
    ibiff_node->value->right = iconn_node;   // CONVENTION: Right will point to subtree
    if (ibiff_node->value->parent != NULL)
    {
        struct segment_node *ibiff_par_node = ibiff_node->value->parent;
        if (ibiff_par_node->value->left->id == iconn_node->id)
            ibiff_par_node->value->left = ibiff_node;
        if (ibiff_par_node->value->right->id == iconn_node->id)
            ibiff_par_node->value->right = ibiff_node;
    }

    // Update ndist from ibiff until the root
    struct segment_node *tmp = ibiff_node;
    while (tmp->value->parent != NULL)
    {
        tmp->value->ndist++;
        tmp = tmp->value->parent;
    }
    tmp->value->ndist++;
    the_network->num_terminals = tmp->value->ndist;

    rescale_tree(ibiff_node,iconn_node,inew_node,Q_perf,delta_p,the_network->num_terminals);

    recalculate_radius(the_network);

    return inew_node;
}

double calc_radius (struct cco_network *the_network, struct segment_node *s)
{
    if (s->value->parent == NULL)
        return s->value->radius;
    else
        return s->value->beta * calc_radius(the_network,s->value->parent);
}

void recalculate_radius (struct cco_network *the_network)
{
    struct segment_list *s_list = the_network->segment_list;
    struct segment_node *tmp = s_list->list_nodes;

    while (tmp != NULL)
    {
        tmp->value->radius = calc_radius(the_network,tmp);

        tmp = tmp->next;
    }
}

void restore_state_tree (struct cco_network *the_network,\
                        struct segment_node *iconn)
{

    double Q_perf = the_network->Q_perf;
    double p_perf = the_network->p_perf;
    double p_term = the_network->p_term;
    double delta_p = p_perf - p_term;

    struct segment_node *ibiff = iconn->value->parent;
    struct segment_node *inew = ibiff->value->left;
    struct segment_node *ibiff_par = ibiff->value->parent;

    // Update "iconn" pointers and values
    iconn->value->parent = ibiff->value->parent;
    iconn->value->src = ibiff->value->src;
    iconn->value->beta = ibiff->value->beta;

    // Update "ibiff" parent subtree pointer if exists
    if (ibiff_par != NULL)
    {
        if (ibiff_par->value->right == ibiff)
            ibiff_par->value->right = iconn;
        else if (ibiff_par->value->left == ibiff)
            ibiff_par->value->left = iconn;
    }   

    // Recalculate R* for "iconn" because we change its length
    if (iconn->value->left == NULL && iconn->value->right == NULL)
        calc_relative_resistance_term(iconn);
    else
        calc_relative_resistance_subtree(iconn,iconn->value->right,iconn->value->left);

    // Eliminate "ibiff" and "inew"
    struct point_list *p_list = the_network->point_list;
    struct segment_list *s_list = the_network->segment_list;

    struct point_node *M = inew->value->src;
    struct point_node *T = inew->value->dest;

    delete_node(p_list,T->id);
    delete_node(p_list,M->id);
    delete_node(s_list,inew->id);
    delete_node(s_list,ibiff->id);

    // Update "ndist" from "iconn" until we reach the root
    struct segment_node *tmp = iconn->value->parent;

    if (tmp != NULL)
    {
        while (tmp->value->parent != NULL)
        {
            tmp->value->ndist--;
            tmp = tmp->value->parent;
        }
        // At the root ...
        tmp->value->ndist--;
        the_network->num_terminals = tmp->value->ndist;
    }
    // Already at the root ...
    else
    {
        the_network->num_terminals = iconn->value->ndist;
    }

    // Update R* and betas until we reach the root
    struct segment_node *ipar = iconn->value->parent;

    // Call the function recursively until we reach the root
    if (ipar != NULL)
    {
        struct segment_node *ipar_left = ipar->value->left;
        struct segment_node *ipar_right = ipar->value->right;

        rescale_until_root(ipar,ipar_left,ipar_right,\
                        Q_perf,delta_p,the_network->num_terminals);
    }
    // Recalculate the root radius when we reach this segment using (2.19)
    else
    {
        iconn->value->radius = pow(iconn->value->resistance * Q_perf / delta_p , 0.25);
    }

    // Update the segments radius
    recalculate_radius(the_network);

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

void calc_pressure_drop_term (struct segment_node *iterm, const double Q_term)
{
    iterm->value->delta_p = iterm->value->resistance * Q_term;
}

void calc_pressure_drop_subtree (struct segment_node *iconn, const double Q_term)
{
    double Q_iconn = iconn->value->ndist * Q_term;

    iconn->value->delta_p = iconn->value->resistance * Q_iconn;
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

double calc_radius_ratio (struct segment_node *iconn, struct segment_node *inew, const double Q_term)
{
    double Q_inew = inew->value->ndist*Q_term;
    double R_inew = inew->value->resistance;

    double Q_iconn = iconn->value->ndist*Q_term;
    double R_iconn = iconn->value->resistance;
    
    double ratio = pow( (Q_iconn * R_iconn) / (Q_inew * R_inew), 0.25 );
    
    return ratio;
}

double calc_bifurcation_ratio (const double radius_ratio, const bool sign)
{
    static const double expoent = -1.0/GAMMA;
    double base = radius_ratio;

    if (sign)
        return pow( 1.0 + ( pow( base, -GAMMA) ) , expoent);
    else
        return pow( 1.0 + ( pow( base, GAMMA) ) , expoent);
}

double calc_segment_volume (struct segment_node *s)
{
    struct point *src = s->value->src->value;
    struct point *dest = s->value->dest->value;

    double l = euclidean_norm(src->x,src->y,src->z,\
                              dest->x,dest->y,dest->z);
    double r = s->value->radius;

    return M_PI * r * r * l;
}

double calc_tree_volume (struct cco_network *the_network)
{
    struct segment_list *s_list = the_network->segment_list;
    struct segment_node *tmp = s_list->list_nodes;
    
    double total_volume = 0.0;

    while (tmp != NULL)
    {
        total_volume += calc_segment_volume(tmp);

        tmp = tmp->next;
    }

    return total_volume;
}

double calc_assymetric_ratio (struct segment_node *right, struct segment_node *left)
{
    // Calculate assymetric ratio using equation (2.12) from Rafael's thesis
    double r_right = right->value->radius;
    double r_left = left->value->radius;

    double epsilon = std::min(r_right,r_left) / std::max(r_right,r_left);
    
    return epsilon;
}

// Output will be given in (ms)
double calc_tau_m (const double cm, const double rm)
{
    // Equation (4.16) from Keener's book 
    return rm * cm;
}

// Output will be given in (cm)
double calc_lambda_m (const double r, const double rc, const double rm)
{
    double d = 2.0 * r;

    // Equation (4.17) from Keener's book
    double num = rm * d;
    double den = 4.0 * rc;

    return sqrt(num / den);
}

// Calculate the propagation velocity using the cable equation (6.21) from Kenner's book
double calc_propagation_velocity (const double r,\
                        const double c, const double cm, const double rc, const double rm)
{
    return c * calc_lambda_m(r,rc,rm) / calc_tau_m(cm,rm);
}

double calc_terminal_activation_time (struct segment_node *s,\
                        const double c, const double cm, const double rc, const double rm)
{
    struct segment_node *tmp = s;
    double at = 0.0;

    while (tmp != NULL)
    {
        at += calc_segment_activation_time(tmp,c,cm,rc,rm);

        tmp = tmp->value->parent;
    }

    return at;
}

double calc_segment_activation_time (struct segment_node *s,\
                        const double c, const double cm, const double rc, const double rm)
{
    struct point *src = s->value->src->value;
    struct point *dest = s->value->dest->value;
    double length = euclidean_norm(src->x,src->y,src->z,\
                                  dest->x,dest->y,dest->z);

    double delta_s = length;
    double r = s->value->radius;

    double velocity = calc_propagation_velocity(r,c,cm,rc,rm);
    //printf("Propagation velocity = %g cm/ms \n",velocity);
    //printf("Distance = %g cm \n",delta_s);

    return delta_s / velocity;
}

// The activation time is already in microseconds
double calc_segment_activation_time_using_level (const double at, struct segment_node *iconn)
{
    double level = calc_segment_level(iconn);

    return pow( at, 1.0/level );
}

double calc_segment_level (struct segment_node *iconn)
{
    struct segment_node *tmp = iconn;

    double level = 1.0;
    while (tmp != NULL)
    {
        level++;

        tmp = tmp->value->parent;
    }
    
    return level;
}

double calc_segment_custom_function_with_level_penalty (const double eval, struct segment_node *iconn)
{
    double level = calc_segment_level(iconn);

    return pow( eval, 1.0/level );
}

double calc_segment_custom_function (struct segment_node *s, const double beta, const double alpha)
{
    struct point *src = s->value->src->value;
    struct point *dest = s->value->dest->value;

    double l = euclidean_norm(src->x,src->y,src->z,\
                              dest->x,dest->y,dest->z);
    double r = s->value->radius;

    return pow(l,beta) * pow(r,alpha);
}

double calc_custom_function (struct cco_network *the_network, const double beta, const double alpha)
{
    struct segment_list *s_list = the_network->segment_list;
    struct segment_node *tmp = s_list->list_nodes;
    
    double total_eval = 0.0;

    while (tmp != NULL)
    {
        total_eval += calc_segment_custom_function(tmp,beta,alpha);

        tmp = tmp->next;
    }

    return total_eval;
}

void make_root_using_cloud_points (struct cco_network *the_network, std::vector<struct point> cloud_points)
{
    int N_term = the_network->N_term;
    double Q_perf = the_network->Q_perf;
    double p_perf = the_network->p_perf;
    double p_term = the_network->p_term;
    double r_perf = the_network->r_perf;
    double A_perf = the_network->A_perf;
    double delta_p = p_perf - p_term;

    int K_term = 1;
    double A_supp = (double)((K_term + 1) * A_perf) / (double)N_term; 
    double r_supp = sqrt(A_supp/M_PI);

    struct point_list *p_list = the_network->point_list;
    struct segment_list *s_list = the_network->segment_list;

    // Positions from the root
    double x_inew[3] = {0,0,0};
    double x_prox[3] = {0,0,0};

    // Sort the distal position of the root until its size is larger than the perfusion radius  
    uint32_t index;
    while (euclidean_norm(x_prox[0],x_prox[1],x_prox[2],x_inew[0],x_inew[1],x_inew[2]) < r_supp)
        index = sort_point_from_cloud(x_inew,cloud_points);
    
    // Eliminate the sorted point from the cloud
    cloud_points.erase(cloud_points.begin()+index);
    
    // Insert points and create the root segment
    struct point_node *A = insert_point(p_list,x_prox);
    struct point_node *B = insert_point(p_list,x_inew);

    struct segment *iroot = new_segment(A,B,NULL,NULL,NULL,Q_perf,p_perf);
    struct segment_node *iroot_node = insert_segment_node(s_list,iroot);
    rescale_root(iroot_node,Q_perf,delta_p);
    the_network->num_terminals = 1;
}

void make_root (struct cco_network *the_network)
{
    int N_term = the_network->N_term;
    double Q_perf = the_network->Q_perf;
    double p_perf = the_network->p_perf;
    double p_term = the_network->p_term;
    double r_perf = the_network->r_perf;
    double A_perf = the_network->A_perf;
    double delta_p = p_perf - p_term;

    int K_term = 1;
    double A_supp = (double)((K_term + 1) * A_perf) / (double)N_term; 
    double r_supp = sqrt(A_supp/M_PI);

    struct point_list *p_list = the_network->point_list;
    struct segment_list *s_list = the_network->segment_list;

    // Positions from the root
    double x_inew[3] = {0,0,0};
    double x_prox[3] = {0,0,0};

    // Sort the distal position of the root until its size is larger than the perfusion radius  
    while (euclidean_norm(x_prox[0],x_prox[1],x_prox[2],x_inew[0],x_inew[1],x_inew[2]) < r_supp)
        generate_point_inside_perfusion_area(x_inew,r_supp);

    // Insert points and create the root segment
    struct point_node *A = insert_point(p_list,x_prox);
    struct point_node *B = insert_point(p_list,x_inew);

    struct segment *iroot = new_segment(A,B,NULL,NULL,NULL,Q_perf,p_perf);
    struct segment_node *iroot_node = insert_segment_node(s_list,iroot);
    rescale_root(iroot_node,Q_perf,delta_p);
    the_network->num_terminals = 1;
}

void generate_terminal_using_cloud_points(struct cco_network *the_network, struct cost_function_config *config,\
                                         std::vector<struct point> cloud_points)
{

    FILE *log_file = the_network->log_file;

    int K_term = the_network->num_terminals;
    int N_term = the_network->N_term;
    double r_perf = the_network->r_perf;
    double A_perf = the_network->A_perf;

    // Cost function reference
    set_cost_function_fn *cost_function_fn = config->function;

    // Increase support domain
    double A_supp = (double)((K_term + 1) * A_perf) / (double)N_term; 
    double r_supp = sqrt(A_supp/M_PI);
    //printf("[!] Support domain radius = %g\n",r_supp);
    fprintf(log_file,"[!] Support domain radius = %g\n",r_supp);

    double new_pos[3];
    bool point_is_ok = false;
    uint32_t tosses = 0;
    double d_threash = calc_dthreashold(r_supp,K_term);

    std::vector<struct segment_node*> feasible_segments;

    // Reference to segment we are going to make the connection
    struct segment_node *iconn = NULL;
    
    while (!point_is_ok)
    {
        // Reset the feasible segments list
        feasible_segments.clear();

        // Sort a terminal position from the cloud of points
        uint32_t index = sort_point_from_cloud(new_pos,cloud_points);

        // RESTRICTION AREA 
        // Check the distance criterion for this point
        point_is_ok = connection_search(the_network,new_pos,d_threash);

        // Check collision with other segments
        if (point_is_ok)
            point_is_ok = check_collisions(the_network,new_pos,feasible_segments);

        // Test the cost function
        if (point_is_ok)
        {
            fprintf(log_file,"Feasible segments: ");
            for (unsigned int i = 0; i < feasible_segments.size(); i++)
                fprintf(log_file,"%d ",feasible_segments[i]->id);
            fprintf(log_file,"\n");

            // COST FUNCTION
            iconn = cost_function_fn(the_network,config,new_pos,feasible_segments);
            if (iconn == NULL)
            {
                fprintf(stderr,"[cco] Error! No feasible segment found!\n");

                point_is_ok = false;
            }
        }

        // If the point does not attend the distance criterion or if there is a collision
        // we need to choose another point.
        if (!point_is_ok)
        {
            tosses++;
            if (tosses > NTOSS)
            {
                //printf("[!] Reducing dthreash! Before = %.2lf || Now = %.2lf \n",\
                        d_threash,d_threash*0.9);
                fprintf(log_file,"[!] Reducing dthreash! Before = %.2lf || Now = %.2lf \n",\
                        d_threash,d_threash*0.9);
                d_threash *= 0.9;
                tosses = 0;
            }
        }
        // The point is valid one and we can eliminate it from the cloud
        else
        {
            cloud_points.erase(cloud_points.begin()+index);
        }
    }

    // Build a new segment with the best connection given by the cost function
    build_segment(the_network,iconn->id,new_pos);
}

void generate_terminal (struct cco_network *the_network, struct cost_function_config *config)
{
    FILE *log_file = the_network->log_file;

    int K_term = the_network->num_terminals;
    int N_term = the_network->N_term;
    double r_perf = the_network->r_perf;
    double A_perf = the_network->A_perf;

    // Cost function reference
    set_cost_function_fn *cost_function_fn = config->function;

    // Increase support domain
    double A_supp = (double)((K_term + 1) * A_perf) / (double)N_term; 
    double r_supp = sqrt(A_supp/M_PI);
    //printf("[!] Support domain radius = %g\n",r_supp);
    fprintf(log_file,"[!] Support domain radius = %g\n",r_supp);

    double new_pos[3];
    bool point_is_ok = false;
    uint32_t tosses = 0;
    double d_threash = calc_dthreashold(r_supp,K_term);
    
    // Array with the reference to the feasible segments for a new connection
    std::vector<struct segment_node*> feasible_segments;

    // Reference to segment we are going to make the connection
    struct segment_node *iconn = NULL;

    while (!point_is_ok)
    {
        // Reset the feasible segments list
        feasible_segments.clear();

        // Generate the terminal position inside the perfusion area
        generate_point_inside_perfusion_area(new_pos,r_supp);

        // RESTRICTION AREA 
        // Check the distance criterion for this point
        point_is_ok = connection_search(the_network,new_pos,d_threash);

        // Check collision with other segments
        if (point_is_ok)
            point_is_ok = check_collisions(the_network,new_pos,feasible_segments);

        // Test the cost function
        if (point_is_ok)
        {
            fprintf(log_file,"Feasible segments: ");
            for (unsigned int i = 0; i < feasible_segments.size(); i++)
                fprintf(log_file,"%d ",feasible_segments[i]->id);
            fprintf(log_file,"\n");

            // COST FUNCTION
            iconn = cost_function_fn(the_network,config,new_pos,feasible_segments);
            if (iconn == NULL)
            {
                //fprintf(stderr,"[cco] Error! No feasible segment found!\n");
                printf("[cco] Error! No feasible segment found!\n");

                point_is_ok = false;
            }
        }

        // If the point does not attend the distance criterion or if there is a collision
        // or if there is no feasible segment for the cost function we need to choose 
        // another point.
        if (!point_is_ok)
        {
            tosses++;
            if (tosses > NTOSS)
            {
                //printf("[!] Reducing dthreash! Before = %.2lf || Now = %.2lf \n",\
                        d_threash,d_threash*0.9);
                fprintf(log_file,"[!] Reducing dthreash! Before = %.2lf || Now = %.2lf \n",\
                        d_threash,d_threash*0.9);
                d_threash *= 0.9;
                tosses = 0;
            }
        }
    }

    // Build a new segment with the best connection given by the cost function
    build_segment(the_network,iconn->id,new_pos);
}

void read_cloud_points (const char filename[], std::vector<struct point> &cloud_points)
{
    printf("[cco] Reading cloud of points !\n");

    uint32_t num_points;
    FILE *file = fopen(filename,"r");
    fscanf(file,"%u",&num_points);

    for (uint32_t i = 0; i < num_points; i++)
    {
        struct point point;
        fscanf(file,"%lf %lf %lf",&point.x,&point.y,&point.z);
    
        cloud_points.push_back(point);
    }

    fclose(file);
}

uint32_t sort_point_from_cloud (double pos[], std::vector<struct point> cloud_points)
{
    uint32_t num_points = cloud_points.size();
    uint32_t index = rand() % num_points;

    pos[0] = cloud_points[index].x;
    pos[1] = cloud_points[index].y;
    pos[2] = cloud_points[index].z;

    return index;
}

void usage (const char pname[])
{
    printf("%s\n",PRINT_LINE);
    printf("Usage:> %s <configuration_file>\n",pname);
    printf("%s\n",PRINT_LINE);
    printf("\t<configuration_file> = Configuration file with the parameter values from the network\n");
    printf("%s\n",PRINT_LINE);
    printf("Example:\n");
    printf("\t%s inputs/simple_cco.ini\n",pname);
    printf("%s\n",PRINT_LINE);
}