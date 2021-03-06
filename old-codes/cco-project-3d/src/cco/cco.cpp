#include "cco.h"

uint32_t cur_rand_index = 0;
uint32_t cur_pmj_index = 0;
double start_radius;

struct cco_network* new_cco_network (struct user_options *options)
{
    struct cco_network *result = (struct cco_network*)malloc(sizeof(struct cco_network));

    set_parameters(result,options);
    set_save_network(result,options);
    set_cost_function_name(result,options);
    set_cloud_points_name(result,options);
    set_obstacle_name(result,options);
    set_pmj_location_name(result,options);
    set_lat_name(result,options);
    set_local_optimization_function_name(result,options);
    set_pruning_function(result,options);

    return result;
}

void free_cco_network (struct cco_network *the_network)
{
    free_point_list(the_network->point_list);
    free_segment_list(the_network->segment_list);
    free(the_network->cost_function_name);
    free(the_network->output_dir);
    if (the_network->using_cloud_points)
        free(the_network->cloud_points_filename);
    if (the_network->using_obstacle)
        free(the_network->obstacle_filename);
    if (the_network->using_local_optimization)
        free(the_network->local_optimization_function_name);
    if (the_network->using_pmj_location)
        free(the_network->pmj_location_filename);
    if (the_network->using_lat)
        free(the_network->lat_filename);
    if (the_network->using_initial_network)
        free(the_network->initial_network_filename);
    if (the_network->using_pruning)
        free(the_network->pruning_function_name);
    fclose(the_network->log_file);

    free(the_network);
}

void set_parameters (struct cco_network *the_network, struct user_options *options)
{
    the_network->point_list = new_point_list();
    the_network->segment_list = new_segment_list();
    the_network->num_terminals = 0;
    the_network->Q_perf = options->q_perf;
    the_network->Q_term = calc_flux_terminals(options->q_perf,options->n_term);
    the_network->p_perf = options->p_perf;
    the_network->p_term = options->p_term;
    the_network->delta_p = options->p_perf - options->p_term;
    the_network->N_term = options->n_term;
    the_network->root_pos[0] = options->root_pos[0];
    the_network->root_pos[1] = options->root_pos[1];
    the_network->root_pos[2] = options->root_pos[2];
    the_network->V_perf = options->v_perf;
    the_network->r_perf = calc_perfusion_radius(options->v_perf);
    the_network->r_supp = the_network->r_perf;                            // When we consider a fixed perfusion volume

    the_network->using_only_murray_law = options->use_only_murray;
    if (options->use_only_murray)
    {
        if (options->start_radius == -1)
        {
            printf("Please enter the initial radius of the root: ");
            scanf("%lf",&start_radius);
        }
        else
        {
            start_radius = options->start_radius;
        }
    }
    the_network->seed = options->seed;
    the_network->max_rand_offset = options->max_rand_offset;
    srand(the_network->seed);

    the_network->using_initial_network = options->use_initial_network;
    if (options->use_initial_network)
    {
        uint32_t size = strlen(options->initial_network_filename) + 1;
        the_network->initial_network_filename = (char*)malloc(sizeof(char)*size);
        strcpy(the_network->initial_network_filename,options->initial_network_filename);
        printf("[cco] Using initial network = %s\n",the_network->initial_network_filename);
    }

    the_network->gamma = options->gamma;
    printf("[cco] Gamma = %g\n",the_network->gamma);
}

void set_cost_function_name (struct cco_network *the_network, struct user_options *options)
{
    uint32_t size = strlen(options->config->function_name) + 1;
    the_network->cost_function_name = (char*)malloc(sizeof(char)*size);
    strcpy(the_network->cost_function_name,options->config->function_name);
    printf("[cco] Cost function name = %s\n",the_network->cost_function_name);
}

void set_save_network (struct cco_network *the_network, struct user_options *options)
{
    uint32_t size = strlen(options->output_dir) + 1;
    the_network->output_dir = (char*)malloc(sizeof(char)*size);
    strcpy(the_network->output_dir,options->output_dir);
    
    create_directory(the_network->output_dir);
    printf("[cco] Output directory = %s\n",the_network->output_dir);

    size += 20;
    char *log_filename = (char*)malloc(sizeof(char)*size);
    sprintf(log_filename,"%s/output.log",the_network->output_dir);
    the_network->log_file = fopen(log_filename,"w+");
    if (!the_network->log_file)
    {
        fprintf(stderr,"[-] ERROR! Opening logfile!\n");
        exit(EXIT_FAILURE);
    }
    free(log_filename);
}

void set_cloud_points_name (struct cco_network *the_network, struct user_options *options)
{
    if (options->use_cloud_points)
    {
        the_network->using_cloud_points = true;
        uint32_t size = strlen(options->cloud_points_filename) + 1;
        the_network->cloud_points_filename = (char*)malloc(sizeof(char)*size);
        strcpy(the_network->cloud_points_filename,options->cloud_points_filename);

        printf("[cco] Using cloud of points\n");
        printf("[cco] Cloud points filename :> \"%s\"\n",the_network->cloud_points_filename);
    }
    else
    {
        the_network->using_cloud_points = false;

        printf("[cco] No cloud of points provided\n");
        printf("[cco] Generating cloud of points based on value of \"r_perf\"\n");
    }
}

void set_obstacle_name (struct cco_network *the_network, struct user_options *options)
{
    if (options->use_obstacle)
    {
        the_network->using_obstacle = true;
        uint32_t size = strlen(options->obstacle_filename) + 1;
        the_network->obstacle_filename = (char*)malloc(sizeof(char)*size);
        strcpy(the_network->obstacle_filename,options->obstacle_filename);

        printf("[cco] Using an obstacle\n");
        printf("[cco] Obstacle points filename :> \"%s\"\n",the_network->obstacle_filename);
    }
    else
    {
        the_network->using_obstacle = false;

        printf("[cco] No obstacle provided\n");
    }
}

void set_pmj_location_name (struct cco_network *the_network, struct user_options *options)
{
    if (options->use_pmj_location)
    {
        the_network->using_pmj_location = true;
        uint32_t size = strlen(options->pmj_location_filename) + 1;
        the_network->pmj_location_filename = (char*)malloc(sizeof(char)*size);
        strcpy(the_network->pmj_location_filename,options->pmj_location_filename);

        printf("[cco] Using PMJ location\n");
        printf("[cco] PMJ location filename :> \"%s\"\n",the_network->pmj_location_filename);
    }
    else
    {
        the_network->using_pmj_location = false;

        printf("[cco] No PMJ location provided\n");
    }
}

void set_lat_name (struct cco_network *the_network, struct user_options *options)
{
    if (options->use_lat)
    {
        the_network->using_lat = true;
        uint32_t size = strlen(options->lat_filename) + 1;
        the_network->lat_filename = (char*)malloc(sizeof(char)*size);
        strcpy(the_network->lat_filename,options->lat_filename);

        printf("[cco] Using LAT from PMJs\n");
        printf("[cco] LAT filename :> \"%s\"\n",the_network->lat_filename);
    }
    else
    {
        the_network->using_lat = false;

        printf("[cco] No LAT provided\n");
    }
}

void set_local_optimization_function_name (struct cco_network *the_network, struct user_options *options)
{
    if (options->use_local_optimization)
    {
        the_network->using_local_optimization = true;
        uint32_t size = strlen(options->local_opt_config->function_name) + 1;
        the_network->local_optimization_function_name = (char*)malloc(sizeof(char)*size);
        strcpy(the_network->local_optimization_function_name,options->local_opt_config->function_name);

        printf("[cco] Using local optimization\n");
        printf("[cco] Local optimization function name :> \"%s\"\n",the_network->local_optimization_function_name);
    }
    else
    {
        the_network->using_local_optimization = false;

        printf("[cco] No local optimization function provided\n");
    }
}

void set_pruning_function (struct cco_network *the_network, struct user_options *options)
{
    the_network->using_pruning = options->use_pruning;
    if (options->use_pruning)
    {
        uint32_t size = strlen(options->pruning_config->function_name) + 1;
        the_network->pruning_function_name = (char*)malloc(sizeof(char)*size);
        strcpy(the_network->pruning_function_name,options->pruning_config->function_name);

        printf("[cco] Using pruning\n");
        printf("[cco] Pruning function name :> \"%s\"\n",the_network->pruning_function_name);
    }
    else
    {
        printf("[cco] No pruning function provided\n");
    }
}

void grow_tree (struct cco_network *the_network, struct user_options *options)
{
    struct stop_watch solver_time;
	init_stop_watch(&solver_time);

    printf("\n[cco] Growing CCO network !\n");

    std::vector<struct point> cloud_points;
    double r_supp = the_network->r_supp;

    // Cloud of points section
    if (!the_network->using_cloud_points)
    {
        build_cloud_points(cloud_points,r_supp);
        draw_perfusion_volume(the_network->r_supp);
    }
    else
        read_cloud_points(the_network->cloud_points_filename,cloud_points);
    
    // Obstacle section
    std::vector<struct face> obstacle_faces;
    if (the_network->using_obstacle)
    {
        read_obstacle_faces(the_network->obstacle_filename,obstacle_faces);
    }

    // PMJ location section
    std::vector<struct point> pmj_points;
    if (the_network->using_pmj_location)
    {
        read_pmj_location_points(the_network->pmj_location_filename,pmj_points);
    }

    // LAT section
    std::vector<struct pmj_lat> lat_points;
    if (the_network->using_lat)
    {
        read_lat_points(the_network->lat_filename,lat_points);
    }
    
    start_stop_watch(&solver_time);

    grow_tree_using_cloud_points(the_network,options,cloud_points,obstacle_faces,pmj_points,lat_points);
    
    long res_time = stop_stop_watch(&solver_time);
	double conv_rate = 1000.0*1000.0*60.0;

    // Unitary test
    check_bifurcation_rule(the_network);
    
    // Print some information about the network
    print_network_info(the_network);
    printf("[INFO] Resolution time = %ld μs (%g min)\n",res_time,res_time/conv_rate);

}

void grow_tree_using_cloud_points (struct cco_network *the_network,\
                                   struct user_options *options,\
                                   std::vector<struct point> cloud_points,\
                                   std::vector<struct face> obstacle_faces,\
                                   std::vector<struct point> pmj_points,\
                                   std::vector<struct pmj_lat> lat_points)
{
    FILE *log_file = the_network->log_file;

    struct cost_function_config *config = options->config;
    struct local_optimization_config *local_opt_config = options->local_opt_config;
    struct pruning_config *pruning_config = options->pruning_config;

    double Q_perf = the_network->Q_perf;
    double p_perf = the_network->p_perf;
    double p_term = the_network->p_term;
    double r_perf = the_network->r_perf;
    double delta_p = p_perf - p_term;

    struct point_list *p_list = the_network->point_list;
    struct segment_list *s_list = the_network->segment_list;

    // ROOT PLACEMENT
    bool using_initial_network = the_network->using_initial_network;
    if (using_initial_network)
    {
        make_root_using_initial_network(the_network);
    }
    else
    {
        make_root_using_cloud_points(the_network,cloud_points,obstacle_faces);
    }   
    
    // MAIN ITERATION LOOP
    while (the_network->num_terminals < the_network->N_term)
    {
        printf("%s\n",PRINT_LINE);
        printf("[cco] Working on terminal number %d\n",the_network->num_terminals+1);
        fprintf(log_file,"%s\n",PRINT_LINE);
        fprintf(log_file,"[cco] Working on terminal number %d\n",the_network->num_terminals);

        generate_terminal_using_cloud_points(the_network,config,local_opt_config,cloud_points,obstacle_faces,lat_points);

        // DEBUG
        write_to_vtk_iteration(the_network);

        printf("%s\n",PRINT_LINE);
        fprintf(log_file,"%s\n",PRINT_LINE);
    }

    // PRUNING
    if (the_network->using_pruning)
    {
        for (uint32_t k = 0; k < PRUNING_PASSES; k++)
            //prune_tree(the_network,pruning_config);
            prune_tree_2(the_network,pruning_config);
    }
        
    // If some terminals were eliminated by the pruning process, add remaining ones
 /*
    while (the_network->num_terminals != the_network->N_term)
    {
        printf("%s\n",PRINT_LINE);
        printf("[cco] Working on terminal number %d after pruning\n",the_network->num_terminals);
        fprintf(log_file,"%s\n",PRINT_LINE);
        fprintf(log_file,"[cco] Working on terminal number %d after pruning\n",the_network->num_terminals);
        printf("%s\n",PRINT_LINE);

        generate_terminal_using_cloud_points(the_network,config,local_opt_config,cloud_points,obstacle_faces,lat_points);
    }
*/
    // PMJ's locations
    if (the_network->using_pmj_location)
    {
        uint32_t num_pmj_points = pmj_points.size();

        for (uint32_t i = 0; i < num_pmj_points; i++)
            generate_terminal_using_pmj_points(the_network,config,local_opt_config,pmj_points,obstacle_faces,lat_points);
    }
        
    // Write to the logfile
    fprintf(log_file,"%s\n",PRINT_LINE);
    write_list(p_list,log_file);
    fprintf(log_file,"%s\n",PRINT_LINE);

    fprintf(log_file,"%s\n",PRINT_LINE);
    write_list(s_list,log_file);
    fprintf(log_file,"%s\n",PRINT_LINE);
}

// p1 -> proximal || p2 -> distal || p3 -> middle || p4 -> new 
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

            bool intersect = collision_detection(src->x,src->y,src->z,\
                                            dest->x,dest->y,dest->z,\
                                            middle_pos[0],middle_pos[1],middle_pos[2],\
                                            new_pos[0],new_pos[1],new_pos[2]);
    
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

bool has_collision (struct segment_list *s_list, struct segment_node *iconn, struct segment_node *ibiff, struct segment_node *inew, FILE *log_file)
{
    // Get the reference to the points from the 'inew' segment
    struct point *src_inew = inew->value->src->value;
    struct point *dest_inew = inew->value->dest->value;

    // Get the reference to the points from the 'ibiff' segment
    struct point *src_ibiff = ibiff->value->src->value;
    struct point *dest_ibiff = ibiff->value->dest->value;

    struct segment_node *tmp = s_list->list_nodes;
    while (tmp != NULL)
    {
        // We avoid both the 'iconn' and 'inew' segments from the search
        if (tmp->id != iconn->id && tmp->id != inew->id && tmp->id != ibiff->id)
        {
            // Get the reference to the points of the current segment
            struct point *src = tmp->value->src->value;
            struct point *dest = tmp->value->dest->value;

            bool intersect = collision_detection(src_inew->x,src_inew->y,src_inew->z,\
                                            dest_inew->x,dest_inew->y,dest_inew->z,\
                                            src->x,src->y,src->z,\
                                            dest->x,dest->y,dest->z);

            if (intersect)
            {
                printf("\t[-] ERROR! Intersection between segments!\n");
                //fprintf(log_file,"\t[-] ERROR! Intersection with segment %d !\n",tmp->id);
                return true;
            }
        }
        tmp = tmp->next;
    }

    return false;
}

bool has_intersect_obstacle (struct segment_node *inew, std::vector<struct face> obstacle_faces)
{
    // Get the reference to the points from the 'inew' segment
    struct point *src_inew = inew->value->src->value;
    struct point *dest_inew = inew->value->dest->value;

    // Get the coordinates from the segment endpoints
    double x_prox[3];
    x_prox[0] = src_inew->x; 
    x_prox[1] = src_inew->y; 
    x_prox[2] = src_inew->z;

    double x_new[3];
    x_new[0] = dest_inew->x; 
    x_new[1] = dest_inew->y; 
    x_new[2] = dest_inew->z;

    return has_intersect_obstacle(x_prox,x_new,obstacle_faces);
}

bool has_intersect_obstacle (const double x_prox[], const double x_new[], std::vector<struct face> obstacle_faces)
{
    for (uint32_t i = 0; i < obstacle_faces.size(); i++)
    {
        if (check_segment_plane_intersection(x_prox,x_new,obstacle_faces[i]))
            return true;
    }

    return false;
}

bool has_valid_segment_sizes (const double iconn_size, const double ibiff_size, const double inew_size)
{
    return (iconn_size > EPSILON && ibiff_size > EPSILON && inew_size > EPSILON);
}

bool has_valid_segment_sizes_2 (const double iconn_size, const double ibiff_size, const double inew_size)
{
    return (iconn_size < MAX_SEGMENT_SIZE && ibiff_size < MAX_SEGMENT_SIZE && inew_size < MAX_SEGMENT_SIZE);
}

bool has_level (struct segment_node *inew, const uint32_t num_terminals)
{
    uint32_t level = calc_segment_level(inew);
    uint32_t level_threashold = calc_level_threashold(num_terminals);

    return (level >= level_threashold);
}

bool check_collisions_and_fill_feasible_segments (struct cco_network *the_network, const double new_pos[],\
                    std::vector<struct segment_node*> &feasible_segments)
{
    FILE *log_file = the_network->log_file;

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

    return (d_crit < d_threash) ? false : true;
}

void rescale_root (struct segment_node *iroot, const double Q_perf, const double delta_p, const bool using_only_murray)
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
    if (!using_only_murray)
    {
        iroot->value->radius = pow(R * Q_perf / delta_p , 0.25);
    }
    // [FRACTAL] The user needs to specify the initial radius of the root
    else
    {
        iroot->value->radius = start_radius;
    }
}

void rescale_tree (struct segment_node *ibiff, struct segment_node *iconn, struct segment_node *inew,\
                 const double Q_perf, const double delta_p, const double gamma, const int num_terminals, const bool use_only_murray)
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
    double radius_ratio;
    if (!use_only_murray)
    {
        radius_ratio = calc_radius_ratio(iconn,inew,Q_term);

        // iconn + inew: Calculate bifurcation ratio using (2.31)
        inew->value->beta = calc_bifurcation_ratio(gamma,radius_ratio,false);
        iconn->value->beta = calc_bifurcation_ratio(gamma,radius_ratio,true);
    }
    // [FRACTAL] Compute the radius ratio using only Murray's law
    else
    {
        // Assymetric
        //double value = 0.5;
        //double r_par = iconn->value->radius;
        //double r_left = value *  powf(0.5,1.0/GAMMA) * r_par;
        //double r_right = powf( powf(r_par,GAMMA) - powf(r_left,GAMMA), 1.0/GAMMA);

        // Symmetric
        double r_par = iconn->value->radius;
        double r_left = pow(0.5, 1.0/gamma) * r_par;
        double r_right = pow(0.5, 1.0/gamma) * r_par;

        //radius_ratio = r_left / r_right;

        // TODO: Use Murray law in some way ...
        // Fix a bifurcation ratio (Decrease the radius at each level of the tree by fixed factor)
        //double level = calc_segment_level(iconn);
        //double ratio;
        //ratio = 0.999;
        inew->value->beta = r_left / r_par;
        iconn->value->beta = r_right / r_par;
    }

    // ibiff: Calculate resistance using (2.5) and pressure drop using (2.7)
    calc_relative_resistance_subtree(ibiff,iconn,inew);
    calc_pressure_drop_subtree(ibiff,Q_term);

    // Rescale the until we reach the root by using the "parent" pointer
    struct segment_node *ipar = ibiff->value->parent;
    if (ipar != NULL)
    {
        struct segment_node *ipar_left = ipar->value->left;
        struct segment_node *ipar_right = ipar->value->right;
        rescale_until_root(ipar,ipar_left,ipar_right,Q_perf,delta_p,gamma,num_terminals,use_only_murray);
    }
    // We are already at the root
    else
    {
        // Recalculate the radius using (2.19)
        if (!use_only_murray)
            ibiff->value->radius = pow(ibiff->value->resistance * Q_perf / delta_p , 0.25);
        // [FRACTAL] Use the initial radius
        else
            ibiff->value->radius = start_radius;
    }
}

void rescale_until_root (struct segment_node *ipar, struct segment_node *ipar_left, struct segment_node *ipar_right,\
                        const double Q_perf, const double delta_p, const double gamma, const int num_terminals, const bool use_only_murray)
{

    // Reach the root
    if (ipar == NULL) return;

    double Q_term = Q_perf / num_terminals;


    if (ipar_left != NULL && ipar_right != NULL)
    {
        //printf("ipar = %u -- ipar_left = %u -- ipar_right = %u\n",ipar->id,ipar_left->id,ipar_right->id);

        double radius_ratio;
        // Calculate radius ratio between left and right subtree
        if (!use_only_murray)
        {
            radius_ratio = calc_radius_ratio(ipar_right,ipar_left,Q_term);

            // Recalculate bifurcation ratios for the offsprings using (2.31)
            ipar_left->value->beta = calc_bifurcation_ratio(gamma,radius_ratio,false);
            ipar_right->value->beta = calc_bifurcation_ratio(gamma,radius_ratio,true);
        }
        // [FRACTAL] Compute the radius ratio using only Murray's law
        else
        {
            // Assymetric
            //double value = 0.5;
            //double r_par = ipar->value->radius;
            //double r_left = value *  powf(0.5,1.0/GAMMA) * r_par;
            //double r_right = powf( powf(r_par,GAMMA) - powf(r_left,GAMMA), 1.0/GAMMA);

            // Symmetric
            double r_par = ipar->value->radius;
            double r_left = pow(0.5, 1.0/gamma) * r_par;
            double r_right = pow(0.5, 1.0/gamma) * r_par;
            ipar_left->value->beta = r_left / r_par;
            ipar_right->value->beta = r_right / r_par;
            
            radius_ratio = r_left / r_right;

            // TODO: Use Murray law in some way ...
            // Fix a bifurcation ratio (Decrease the radius at each level of the tree by fixed factor)
            //double level = calc_segment_level(ipar_right);
            //double ratio;
            //ratio = 0.999;
            //ipar_left->value->beta = ratio;
            //ipar_right->value->beta = ratio;
        }

        // Recalculate resistance using (2.5) and pressure drop using (2.7)
        calc_relative_resistance_subtree(ipar,ipar_left,ipar_right);
        calc_pressure_drop_subtree(ipar,Q_term);

        // Call the function recursively until we reach the root
        if (ipar->value->parent != NULL)
            rescale_until_root(ipar->value->parent,\
                               ipar->value->parent->value->left,\
                               ipar->value->parent->value->right,\
                               Q_perf,delta_p,gamma,num_terminals,use_only_murray);
        else
        {
            // Recalculate the root radius when we reach this segment using (2.19)
            if (!use_only_murray)
                ipar->value->radius = pow(ipar->value->resistance * Q_perf / delta_p , 0.25);
            // [FRACTAL] Use the initial root radius
            else
                ipar->value->radius = start_radius;
        }
    }
}

struct segment_node* build_segment (struct cco_network *the_network, struct local_optimization_config *local_opt_config,\
                                const uint32_t index, const double new_pos[])
{

    struct point_list *p_list = the_network->point_list;
    struct segment_list *s_list = the_network->segment_list;
    double Q_perf = the_network->Q_perf;
    double p_perf = the_network->p_perf;
    double p_term = the_network->p_term;
    double delta_p = p_perf - p_term;
    double gamma = the_network->gamma;
    bool using_local_optimization = the_network->using_local_optimization;

    // TODO: Pass the pointer directly here ...
    struct segment_node *iconn_node = search_segment_node(s_list,index);
    struct segment *iconn = iconn_node->value;

    // Create the bifurcation point based on the user settings
    struct point_node *M;
    if (using_local_optimization)
    {
        bool first_call = local_opt_config->first_call;
    
        if (first_call)
        {
            double middle_pos[3];
            calc_middle_point_segment(iconn_node,middle_pos);  
            M = insert_point(p_list,middle_pos); 
        }
        else
        {
            double *best_pos = local_opt_config->best_pos;
            M = insert_point(p_list,best_pos);

            // Set on the 'first_call' flag
            local_opt_config->first_call = true;
        }
        
    }
    else
    {
        double middle_pos[3];
        calc_middle_point_segment(iconn_node,middle_pos);
        M = insert_point(p_list,middle_pos);
    }

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

    rescale_tree(ibiff_node,iconn_node,inew_node,Q_perf,delta_p,gamma,the_network->num_terminals,the_network->using_only_murray_law);

    recalculate_radius(the_network);

    return inew_node;
}

double calc_radius (struct cco_network *the_network, struct segment_node *s)
{
    //printf("beta = %g\n",s->value->beta);
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

    double gamma = the_network->gamma;
    double Q_perf = the_network->Q_perf;
    double p_perf = the_network->p_perf;
    double p_term = the_network->p_term;
    double delta_p = p_perf - p_term;

    bool use_only_murray = the_network->using_only_murray_law;

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
                        Q_perf,delta_p,gamma,the_network->num_terminals,use_only_murray);
        
    }
    else
    {
        // Recalculate the root radius when we reach this segment using (2.19)
        if (!use_only_murray)
        {
            iconn->value->radius = pow(iconn->value->resistance * Q_perf / delta_p , 0.25);
        }
        // [FRACTAL] Use the initial value of the root radius
        else
        {
            iconn->value->radius = start_radius;
        }
    }

    // Update the segments radius
    recalculate_radius(the_network);

}

void make_root_using_cloud_points (struct cco_network *the_network, std::vector<struct point> cloud_points, std::vector<struct face> obstacle_faces)
{
    int N_term = the_network->N_term;
    double Q_perf = the_network->Q_perf;
    double p_perf = the_network->p_perf;
    double p_term = the_network->p_term;
    double r_perf = the_network->r_perf;
    double V_perf = the_network->V_perf;
    double delta_p = p_perf - p_term;
    double *root_pos = the_network->root_pos;
    bool using_only_murray_law = the_network->using_only_murray_law;

    uint32_t max_rand_offset = the_network->max_rand_offset;

    struct point_list *p_list = the_network->point_list;
    struct segment_list *s_list = the_network->segment_list;

    int K_term = 1;
    double r_supp = r_perf;

    // Set the position of the root
    double x_prox[3], x_inew[3];
    for (uint32_t i = 0; i < 3; i++)
    {
        x_prox[i] = root_pos[i];
        x_inew[i] = 0.0;
    }

    // Create the first node of tree with proximal root position
    struct point_node *A = insert_point(p_list,x_prox);

    // Sort the distal position of the root using the cloud of points until the root segment has a size larger than a threashold distance
    uint32_t index = 0;
    uint32_t counter = 0;
    uint32_t iterations = 0;
    bool is_root_ok = false;
    double d_threashold = sqrt( M_PI * r_supp * r_supp );
    while (!is_root_ok)
    {
        //sort_point_from_cloud_v1(x_inew,cloud_points);
        //sort_point_from_cloud_v2(x_inew,cloud_points);
        //sort_point_from_cloud_v3(x_inew,cloud_points);
        sort_point_from_cloud_v4(x_inew,cloud_points,max_rand_offset);

        // Convert to the real domain
        //x_inew[0] *= r_supp;
        //x_inew[1] *= r_supp;
        //x_inew[2] *= r_supp;
        double root_length = euclidean_norm(x_prox[0],x_prox[1],x_prox[2],x_inew[0],x_inew[1],x_inew[2]);

        if (root_length >= d_threashold && !has_intersect_obstacle(x_prox,x_inew,obstacle_faces))
            is_root_ok = true;
        else
            counter++;
        
        if (counter > 8)
            d_threashold *= 0.9;

        iterations++;
    }
    printf("[cco] Root segment was set in %u iterations.\n", iterations);

    // Create the distal point of the root
    struct point_node *B = insert_point(p_list,x_inew);

    // Create the root segment
    struct segment *iroot = new_segment(A,B,NULL,NULL,NULL,Q_perf,p_perf);
    struct segment_node *iroot_node = insert_segment_node(s_list,iroot);
    rescale_root(iroot_node,Q_perf,delta_p,using_only_murray_law);
    the_network->num_terminals = 1;

}

void make_root_using_initial_network (struct cco_network *the_network)
{
    printf("[cco] Building root using initial network '%s' ...\n",the_network->initial_network_filename);

    int N_term = the_network->N_term;
    double Q_perf = the_network->Q_perf;
    double p_perf = the_network->p_perf;
    double p_term = the_network->p_term;
    double r_perf = the_network->r_perf;
    double V_perf = the_network->V_perf;
    double delta_p = p_perf - p_term;

    struct point_list *p_list = the_network->point_list;
    struct segment_list *s_list = the_network->segment_list;    

    std::vector<struct cco_point> the_points;
    std::vector<struct cco_segment> the_segments;
    read_initial_network(the_network->initial_network_filename,the_points,the_segments);

    // Insert the points
    for (uint32_t i = 0; i < the_points.size(); i++)
    {
        double pos[3];
        pos[0] = the_points[i].x;
        pos[1] = the_points[i].y;
        pos[2] = the_points[i].z;

        struct point_node *new_p = insert_point(p_list,pos);
    }
    // Insert the segments
    for (uint32_t i = 0; i < the_segments.size(); i++)
    {
        int ids[2];
        ids[0] = the_segments[i].src_id;
        ids[1] = the_segments[i].dest_id;

        struct point_node *src = search_node(p_list,ids[0]);
        struct point_node *dest = search_node(p_list,ids[1]);

        struct segment *new_seg = new_segment(src,dest,NULL,NULL,NULL,Q_perf,p_perf);
        struct segment_node *new_seg_node = insert_segment_node(s_list,new_seg);
    }
    // Set segment pointers
    struct segment_node *tmp = s_list->list_nodes;
    for (uint32_t i = 0; i < the_segments.size(); i++)
    {
        int links[3];
        links[0] = the_segments[i].parent_id;
        links[1] = the_segments[i].left_offspring_id;
        links[2] = the_segments[i].right_offspring_id;

        struct segment_node *parent = NULL;
        struct segment_node *left = NULL;
        struct segment_node *right = NULL;

        if (links[0] != -1)
            parent = search_segment_node(s_list,links[0]);
        if (links[1] != -1)
            left = search_segment_node(s_list,links[1]);
        if (links[2] != -1)
            right = search_segment_node(s_list,links[2]);

        tmp->value->parent = parent;
        tmp->value->left = left;
        tmp->value->right = right;

        tmp = tmp->next;
    }
    // Update "ndist" from every terminal
    tmp = s_list->list_nodes;
    for (uint32_t i = 0; i < the_segments.size(); i++)
    {
        tmp->value->ndist = 0;
        tmp->value->radius = 0.00102269;            // TODO: This information must come from the input file
        if (is_terminal(tmp))
        {
            struct segment_node *tmp_2 = tmp;

            if (tmp_2 != NULL)
            {
                while (tmp_2->value->parent != NULL)
                {
                    tmp_2->value->ndist++;
                    tmp_2 = tmp_2->value->parent;
                }
                // At the root ...
                tmp_2->value->ndist++;
            }
            the_network->num_terminals++;
        }
        tmp = tmp->next;
    }
}

void generate_terminal_using_cloud_points(struct cco_network *the_network,\
                                          struct cost_function_config *config,\
                                          struct local_optimization_config *local_opt_config,\
                                          std::vector<struct point> cloud_points,\
                                          std::vector<struct face> obstacle_faces,\
                                          std::vector<struct pmj_lat> lat_points)
{

    FILE *log_file = the_network->log_file;

    bool using_local_optimization = the_network->using_local_optimization;
    int K_term = the_network->num_terminals;
    int N_term = the_network->N_term;
    double r_perf = the_network->r_perf;
    double r_supp = the_network->r_supp;
    double V_perf = the_network->V_perf;

    uint32_t max_rand_offset = the_network->max_rand_offset;

    // Cost function reference
    set_cost_function_fn *cost_function_fn = config->function;

    // Local optimization reference
    set_local_optimization_function_fn *local_optimization_fn = local_opt_config->function;

    bool is_point_ok = false;
    bool connect_point = true;

    // Calculate the distance threashold
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
        //sort_point_from_cloud_v1(new_pos,cloud_points);
        //sort_point_from_cloud_v2(new_pos,cloud_points);
        //sort_point_from_cloud_v3(new_pos,cloud_points);
        sort_point_from_cloud_v4(new_pos,cloud_points,max_rand_offset);

        // Convert to the real domain
        //new_pos[0] *= r_supp;
        //new_pos[1] *= r_supp;
        //new_pos[2] *= r_supp;

        // RESTRICTION AREA
        // Check the distance criterion for this point
        point_is_ok = connection_search(the_network,new_pos,d_threash);
        //printf("%d -- dist = %d -- (%g,%g,%g)\n",cur_rand_index,point_is_ok,new_pos[0],new_pos[1],new_pos[2]);

        // Check collision with other segments
        if (point_is_ok)
            point_is_ok = check_collisions_and_fill_feasible_segments(the_network,new_pos,feasible_segments);

        // Test the cost function
        if (point_is_ok)
        {
            fprintf(log_file,"Feasible segments: ");
            for (unsigned int i = 0; i < feasible_segments.size(); i++)
                fprintf(log_file,"%d ",feasible_segments[i]->id);
            fprintf(log_file,"\n");

            //printf("total bifurcations = %u\n",feasible_segments.size());

            // COST FUNCTION
            iconn = cost_function_fn(the_network,config,local_opt_config,new_pos,feasible_segments,obstacle_faces,lat_points);
            if (iconn == NULL)
            {
                //fprintf(stderr,"[cco] Error! No feasible segment found!\n");
                //exit(EXIT_FAILURE);

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
                fprintf(log_file,"[!] Reducing dthreash! Before = %g || Now = %g \n",\
                        d_threash,d_threash*FACTOR);
                d_threash *= FACTOR;
                tosses = 0;
            }
        }
        // The point is a valid one and we can eliminate it from the cloud
        else
        {
            //cloud_points.erase(cloud_points.begin()+index);
        }
    }

    // Build the new segment of the tree
    build_segment(the_network,local_opt_config,iconn->id,new_pos);
}

void generate_terminal_using_pmj_points(struct cco_network *the_network,\
                                          struct cost_function_config *config,\
                                          struct local_optimization_config *local_opt_config,\
                                          std::vector<struct point> pmj_points,\
                                          std::vector<struct face> obstacle_faces,\
                                          std::vector<struct pmj_lat> lat_points)
{

    FILE *log_file = the_network->log_file;

    bool using_local_optimization = the_network->using_local_optimization;
    int K_term = the_network->num_terminals;
    int N_term = the_network->N_term;
    double r_perf = the_network->r_perf;
    double r_supp = the_network->r_supp;
    double V_perf = the_network->V_perf;

    // Cost function reference
    set_cost_function_fn *cost_function_fn = config->function;

    // Local optimization reference
    set_local_optimization_function_fn *local_optimization_fn = local_opt_config->function;

    bool is_point_ok = false;
    bool connect_point = true;

    // Calculate the distance threashold
    double new_pos[3];
    bool point_is_ok = false;
    uint32_t tosses = 0;
    double d_threash = calc_dthreashold(r_supp,K_term);
    
    std::vector<struct segment_node*> feasible_segments;

    // Reference to segment we are going to make the connection
    struct segment_node *iconn = NULL;

    // Get a PMJ position
    printf("[cco] Trying to insert PMJ point %u ...\n",cur_pmj_index+1);
    fprintf(log_file,"[cco] Trying to insert PMJ point %u ...\n",cur_pmj_index+1);
    new_pos[0] = pmj_points[cur_pmj_index].x;
    new_pos[1] = pmj_points[cur_pmj_index].y;
    new_pos[2] = pmj_points[cur_pmj_index].z;
    cur_pmj_index++;

    while (!point_is_ok)
    {
        // Reset the feasible segments list
        feasible_segments.clear();

        // Convert to the real domain
        //new_pos[0] *= r_supp;
        //new_pos[1] *= r_supp;
        //new_pos[2] *= r_supp;

        // RESTRICTION AREA
        // Check the distance criterion for this point
        point_is_ok = connection_search(the_network,new_pos,d_threash);

        // Check collision with other segments
        if (point_is_ok)
            point_is_ok = check_collisions_and_fill_feasible_segments(the_network,new_pos,feasible_segments);

        // Test the cost function
        if (point_is_ok)
        {
            fprintf(log_file,"Feasible segments: ");
            for (unsigned int i = 0; i < feasible_segments.size(); i++)
                fprintf(log_file,"%d ",feasible_segments[i]->id);
            fprintf(log_file,"\n");

            //printf("total bifurcations = %u\n",feasible_segments.size());

            // COST FUNCTION
            iconn = cost_function_fn(the_network,config,local_opt_config,new_pos,feasible_segments,obstacle_faces,lat_points);
            if (iconn == NULL)
            {
                //fprintf(stderr,"[cco] Error! No feasible segment found!\n");
                //exit(EXIT_FAILURE);

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
                fprintf(log_file,"[!] Reducing dthreash! Before = %g || Now = %g \n",\
                        d_threash,d_threash*FACTOR);
                d_threash *= FACTOR;
                tosses = 0;
            }

            // GIVE UP !!!!
            if (d_threash < 1.0e-100)
                return;
        }
        // The point is a valid one and we can eliminate it from the cloud
        else
        {
            //cloud_points.erase(cloud_points.begin()+index);
        }
    }

    // Build the new segment of the tree
    build_segment(the_network,local_opt_config,iconn->id,new_pos);
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

void read_pmj_location_points (const char filename[], std::vector<struct point> &pmj_points)
{
    printf("[cco] Reading PMJ locations from '%s' !\n",filename);

    char str[200];
    FILE *file = fopen(filename,"r");

    uint32_t num_points;
    fscanf(file,"%u",&num_points);

    for (uint32_t i = 0; i < num_points; i++)
    {
        struct point point;
        fscanf(file,"%lf %lf %lf",&point.x,&point.y,&point.z);

        pmj_points.push_back(point);
    }

    fclose(file);
}

void read_lat_points (const char filename[], std::vector<struct pmj_lat> &lat_points)
{
    printf("[cco] Reading LAT from '%s' !\n",filename);

    char str[200];
    FILE *file = fopen(filename,"r");

    uint32_t num_points;
    fscanf(file,"%u",&num_points);

    for (uint32_t i = 0; i < num_points; i++)
    {
        struct pmj_lat point;
        fscanf(file,"%lf %lf %lf %lf",&point.x,&point.y,&point.z,&point.at);
        point.is_active = true;

        lat_points.push_back(point);
    }

    fclose(file);
}

void read_obstacle_faces (const char filename[], std::vector<struct face> &obstacle_faces)
{
    printf("[cco] Reading obstacle faces from '%s' !\n",filename);

    char str[200];
    FILE *file = fopen(filename,"r");

    while (fscanf(file,"%s",str) != EOF)
    {
        if (strcmp(str,"facet") == 0)
            read_face(file,obstacle_faces);
    }
    fclose(file);
}

void read_face (FILE *file, std::vector<struct face> &faces)
{
    char str[200];
    double n[3];
    double a[3], b[3], c[3];

    // Read normal vector
    fscanf(file,"%s",str);
    fscanf(file,"%lf %lf %lf",&n[0],&n[1],&n[2]);

    // Read vertex
    fscanf(file,"%s",str); fscanf(file,"%s",str);
    fscanf(file,"%s",str);
    fscanf(file,"%lf %lf %lf",&a[0],&a[1],&a[2]);
    //printf("%g %g %g\n",a[0],a[1],a[2]);
    fscanf(file,"%s",str);
    fscanf(file,"%lf %lf %lf",&b[0],&b[1],&b[2]);
    //printf("%g %g %g\n",b[0],b[1],b[2]);
    fscanf(file,"%s",str);
    fscanf(file,"%lf %lf %lf",&c[0],&c[1],&c[2]);
    //printf("%g %g %g\n\n",c[0],c[1],c[2]);
    fscanf(file,"%s",str); fscanf(file,"%s",str);

    struct face f;
    f.x1 = a[0]; f.y1 = a[1]; f.z1 = a[2];
    f.x2 = b[0]; f.y2 = b[1]; f.z2 = b[2];
    f.x3 = c[0]; f.y3 = c[1]; f.z3 = c[2];
    f.nx = n[0]; f.ny = n[1]; f.nz = n[2];
    
    faces.push_back(f);
}

void build_cloud_points (std::vector<struct point> &cloud_points, const double radius)
{
    printf("[cco] Building cloud of points !\n");

    struct random_generator *the_generator = new_random_generator();
    generate_random_array(the_generator);

    generate_cloud_points(the_generator,cloud_points,radius);

    free_random_generator(the_generator);
}

// Randomly choose a point inside the cloud
void sort_point_from_cloud_v1 (double pos[], std::vector<struct point> cloud_points)
{
    uint32_t num_points = cloud_points.size();
    uint32_t index = rand() % num_points;

    pos[0] = cloud_points[index].x;
    pos[1] = cloud_points[index].y;
    pos[2] = cloud_points[index].z;

    //return index;
}

// Sort points in a sequential order
void sort_point_from_cloud_v2 (double pos[], std::vector<struct point> cloud_points)
{
    // Reset the counter
    if (cur_rand_index == cloud_points.size()-1)
        cur_rand_index = 0;

    // Convert to the real domain
    pos[0] = cloud_points[cur_rand_index].x;
    pos[1] = cloud_points[cur_rand_index].y;
    pos[2] = cloud_points[cur_rand_index].z;

    // Increase the counter
    cur_rand_index++;
}

// Sort points in a sequential order, but using a fixed offset
void sort_point_from_cloud_v3 (double pos[], std::vector<struct point> cloud_points)
{
    static const uint32_t offset = 1;     
    
    // Reset the counter
    if (cur_rand_index > cloud_points.size()-1)
        cur_rand_index = cur_rand_index % (cloud_points.size()-1);

    // Convert to the real domain
    pos[0] = cloud_points[cur_rand_index].x;
    pos[1] = cloud_points[cur_rand_index].y;
    pos[2] = cloud_points[cur_rand_index].z;

    // Increase the counter
    cur_rand_index += offset;
}

// Sort points in a sequential order, but using a variable offset
void sort_point_from_cloud_v4 (double pos[], std::vector<struct point> cloud_points, uint32_t max_rand_offset)
{
    uint32_t offset = rand() % max_rand_offset + 1;           // max_offset = 5 --> Nice results !!!!

    // Reset the counter
    if (cur_rand_index > cloud_points.size()-1)
        cur_rand_index = cur_rand_index % (cloud_points.size()-1);

    // Convert to the real domain
    pos[0] = cloud_points[cur_rand_index].x;
    pos[1] = cloud_points[cur_rand_index].y;
    pos[2] = cloud_points[cur_rand_index].z;

    // Increase the counter
    cur_rand_index += offset;
}

void prune_tree_segment (struct cco_network *the_network, struct segment_node *inew)
{

    double gamma = the_network->gamma;
    double Q_perf = the_network->Q_perf;
    double p_perf = the_network->p_perf;
    double p_term = the_network->p_term;
    double delta_p = p_perf - p_term;

    bool use_only_murray = the_network->using_only_murray_law;

    struct segment_node *ibiff;
    struct segment_node *ibiff_par;
    struct segment_node *iconn;
    
    ibiff = inew->value->parent;
    if (ibiff->value->left == inew)
        iconn = ibiff->value->right;
    else if (ibiff->value->right == inew)
        iconn = ibiff->value->left;
    ibiff_par = ibiff->value->parent;

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
                        Q_perf,delta_p,gamma,the_network->num_terminals,use_only_murray);
    }
    else
    {
        // Recalculate the root radius when we reach this segment using (2.19)
        if (!use_only_murray)
        {
            iconn->value->radius = pow(iconn->value->resistance * Q_perf / delta_p , 0.25);
        }
        // [FRACTAL] Use the initial value of the root radius
        else
        {
            iconn->value->radius = start_radius;
        }
    }

    // Update the segments radius
    recalculate_radius(the_network);

}

void prune_tree (struct cco_network *the_network, struct pruning_config *config)
{
    // User parameters for the pruning function
    set_pruning_function_fn *pruning_fn = config->function;
    
    struct segment_node *tmp = the_network->segment_list->list_nodes;
    while (tmp != NULL)
    {
        // Save the reference to the next node
        struct segment_node *tmp_next = tmp->next;

        // We only prune terminal segments
        if (is_terminal(tmp) && tmp->id < the_network->segment_list->num_nodes)
        {
            // Calculate the segment level and length
            double level = calc_segment_level(tmp);
            double length = calc_segment_size(tmp);
            
            // Call the pruning function
            double eval = pruning_fn(config,level,length);
            
            double r = ((double)rand() / (double)RAND_MAX) * 100.0;

            if (r <= eval)
            {
                printf("[!] Pruning segment %u!\n",tmp->id);
                prune_tree_segment(the_network,tmp);
            }       
        }

        tmp = tmp_next;
    }
}

void prune_tree_2 (struct cco_network *the_network, struct pruning_config *config)
{
    // User parameters for the pruning function
    set_pruning_function_fn *pruning_fn = config->function;
    
    struct segment_node *tmp = the_network->segment_list->list_nodes;
    while (tmp != NULL)
    {

        // We only prune terminal segments
        if (is_terminal(tmp) && tmp->id < the_network->segment_list->num_nodes)
        {
            // Calculate the segment level and length
            double level = calc_segment_level(tmp);
            double length = calc_segment_size(tmp);
            
            // Call the pruning function
            double eval = pruning_fn(config,level,length);
            
            double r = ((double)rand() / (double)RAND_MAX) * 100.0;

            if (r <= eval)
            {
                printf("[!] Pruning segment %u!\n",tmp->id);
                tmp->value->prune = true;
            }       
        }
        tmp = tmp->next;
    }

    tmp = the_network->segment_list->list_nodes;
    while (tmp != NULL)
    {
        if (tmp->value->prune)
        {
            printf("[!] Segment %u has been pruned!\n",tmp->id);
            prune_tree_segment(the_network,tmp);
            tmp = the_network->segment_list->list_nodes;
        }
        tmp = tmp->next;
    }
}

void print_network_info (struct cco_network *the_network)
{
    struct segment_node *tmp = the_network->segment_list->list_nodes;
    double mean_biff_angle = 0.0;
    double mean_segment_length = 0.0;
    uint32_t biff_counter = 0;
    uint32_t num_segments = the_network->segment_list->num_nodes;

    size_t size = strlen(the_network->output_dir) + 100;
    char *filename = (char*)malloc(sizeof(char)*size);

    FILE *file_segments_length, *file_bifurcations;
    double u[3], v[3];
    double angle;

    sprintf(filename,"%s/segments_length.dat",the_network->output_dir);
    file_segments_length = fopen(filename,"w+");
    sprintf(filename,"%s/bifurcations_angle.dat",the_network->output_dir);
    file_bifurcations = fopen(filename,"w+");

    while (tmp != NULL)
    {
        // Segment length calculus
        double length = euclidean_norm(tmp->value->src->value->x,tmp->value->src->value->y,tmp->value->src->value->z,\
                                        tmp->value->dest->value->x,tmp->value->dest->value->y,tmp->value->dest->value->z);
        mean_segment_length += length;
        fprintf(file_segments_length,"%g\n",length);
        
        if (tmp->value->left != NULL && tmp->value->right != NULL)
        {   
            // Bifurcation angle calculus
            calc_unitary_vector(tmp->value->left,u);
            calc_unitary_vector(tmp->value->right,v);
            angle = calc_angle_between_vectors(u,v);
            fprintf(file_bifurcations,"%g\n",angle);

            mean_biff_angle += angle;
            biff_counter++;
        }
        tmp = tmp->next;
    }
    mean_segment_length /= (double)num_segments;
    mean_biff_angle /= (double)biff_counter;

    // Convert the length to {mm}
    mean_segment_length *= 1000;

    fclose(file_segments_length);
    fclose(file_bifurcations);

    printf("[INFO] Total number of segment = %u\n",num_segments);
    printf("[INFO] Mean segment length = %g mm\n",mean_segment_length);
    printf("[INFO] Total number of bifurcations = %u\n",biff_counter);
    printf("[INFO] Mean bifurcation angle = %g degrees\n",mean_biff_angle);

    sprintf(filename,"%s/network_info.dat",the_network->output_dir);
    FILE *file_info = fopen(filename,"w+");

    fprintf(file_info,"total_num_segments %u\n",num_segments);
    fprintf(file_info,"mean_segment_length %g\n",mean_segment_length);
    fprintf(file_info,"total_num_bifurcations %u\n",biff_counter);
    fprintf(file_info,"mean_bifurcation_angle %g\n",mean_biff_angle);

    free(filename);
    fclose(file_info);

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
