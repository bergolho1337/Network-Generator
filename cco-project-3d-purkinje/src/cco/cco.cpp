#include "cco.h"

// Constant definitions
const double CCO_Network::Q_PERF = 8.33e-06;
const double CCO_Network::P_PERF = 13300;
const double CCO_Network::P_TERM = 9580;
const double CCO_Network::V_PERF = 0.0001;

CCO_Network::CCO_Network () { }

CCO_Network::CCO_Network (User_Options *options)
{
    set_parameters(options);
    set_save_network();
    set_cost_function();
    set_local_optimization_function_name();
    
    //print();
}

CCO_Network::~CCO_Network ()
{
    if (this->cost_fn)
        delete this->cost_fn;
    if (this->local_opt_fn)
        delete this->local_opt_fn;
    if (this->pmj_data)
        delete this->pmj_data;
    if (this->cloud_data)
        delete this->cloud_data;
    for (uint32_t i = 0; i < this->segment_list.size(); i++)
        delete this->segment_list[i];
    for (uint32_t i = 0; i < this->point_list.size(); i++)
        delete this->point_list[i];
    if (this->log_file)
        fclose(this->log_file);
}

void CCO_Network::grow_tree (User_Options *options)
{
    struct stop_watch solver_time;
	init_stop_watch(&solver_time);

    printf("\n[cco] Growing CCO network !\n");

    this->cloud_data = new Cloud(options->cloud_config);
    this->pmj_data = new PMJ(options->pmj_config);
    //read_obstacles();

    start_stop_watch(&solver_time);

    grow_tree_using_cloud_points(options);

    long res_time = stop_stop_watch(&solver_time);
	double conv_rate = 1000.0*1000.0*60.0;
    printf("[INFO] Resolution time = %ld μs (%g min)\n",res_time,res_time/conv_rate);

    // Unitary test
    //bool sucess = check_bifurcation_rule(this->log_file,this->gamma,this->segment_list);
    //if (!sucess)
    //{
    //    printf("[-] ERROR! Bifurcation rule failure!\n");
    //    exit(EXIT_FAILURE);
    //}

    // Output network information
    print_network_info();
}

void CCO_Network::set_parameters (User_Options *options)
{
    this->num_terminals = 0;
    this->N_term = options->n_term;
    this->seed = options->seed;
    this->max_rand_offset = options->max_rand_offset;
    this->gamma = options->gamma;
    this->lat_offset = options->lat_offset;
    this->cur_pmj_package = 0;
    this->cur_rand_index = 0;
    this->max_lat_error = __DBL_MIN__;
    this->min_max_aprox_lat[0] = __DBL_MAX__; this->min_max_aprox_lat[1] = __DBL_MIN__;
    this->min_max_ref_lat[0] = __DBL_MAX__; this->min_max_ref_lat[1] = __DBL_MIN__;
    this->min_term_lat = __DBL_MAX__; this->max_term_lat = __DBL_MIN__;
    srand(this->seed);

    memcpy(this->root_pos,options->root_pos,sizeof(double)*3);

    this->using_only_murray_law = options->use_only_murray;
    if (options->use_only_murray)
    {
        if (options->start_radius == -1)
        {
            printf("Please enter the initial radius of the root: ");
            scanf("%lf",this->start_radius);
        }
        else
        {
            this->start_radius = options->start_radius;
        }
    }

    this->using_cloud_points = options->use_cloud_points;
    this->using_local_optimization = options->use_local_optimization;
    this->using_pmj_location = options->use_pmj_location;

    this->output_dir = options->output_dir;
    this->cost_function_name = options->cost_function_config->function_name;
    this->local_optimization_function_name = options->local_opt_config->function_name;

    this->Q_TERM = Q_PERF / this->N_term;
    this->R_PERF = powf(((3.0*V_PERF)/(4.0 * M_PI)), (1.0/3.0));
}

void CCO_Network::set_cost_function ()
{
    if (this->cost_function_name == "minimize_custom_function")
    {
        this->cost_fn = new MinimizeCustomFunction;
    }
    else if (this->cost_function_name == "maximize_custom_function")
    {
        this->cost_fn = new MaximizeCustomFunction;
    }
    else
    {
        fprintf(stderr,"[-] ERROR! Cost function '%s' not found!\n",this->cost_function_name.c_str());
        fprintf(stderr,"Available options are:\n");
        fprintf(stderr,"\t1) minimize_custom_function\n");
        fprintf(stderr,"\t2) minimize_activation_time_function\n");
        exit(EXIT_FAILURE);
    }
}

void CCO_Network::set_local_optimization_function_name ()
{
    if (this->local_optimization_function_name == "default_local_optimization")
    {
        this->local_opt_fn = new DefaultOptimization;
    }
    else if (this->local_optimization_function_name == "rafael_local_optimization")
    {
        this->local_opt_fn = new RafaelOptimization;
    }
    else
    {
        fprintf(stderr,"[-] ERROR! Local optimization function '%s' not found!\n",this->local_optimization_function_name.c_str());
        fprintf(stderr,"Available options are:\n");
        fprintf(stderr,"\t1) default_local_optimization\n");
        fprintf(stderr,"\t2) rafael_local_optimization\n");
        exit(EXIT_FAILURE);
    }
}

void CCO_Network::set_save_network ()
{
    create_directory(this->output_dir.c_str());
    printf("[cco] Output directory = %s\n",this->output_dir.c_str());

    std::string log_filename = this->output_dir + "/output.log";
    this->log_file = fopen(log_filename.c_str(),"w+");
    if (!this->log_file)
    {
        fprintf(stderr,"[-] ERROR! Opening logfile!\n");
        exit(EXIT_FAILURE);
    }
}

void CCO_Network::grow_tree_using_cloud_points (User_Options *options)
{
    bool sucess;

    CostFunctionConfig *cost_function_config = options->cost_function_config;
    LocalOptimizationConfig *local_opt_config = options->local_opt_config;
    //struct pruning_config *pruning_config = options->pruning_config;

    // Root placement
    bool using_initial_network = options->use_initial_network;
    (using_initial_network) ? make_root_using_initial_network() : make_root_using_cloud_points();
    
    // ===============================================================================================================
    // MAIN ITERATION LOOP
    while (this->num_terminals < this->N_term)
    {
        printf("%s\n",PRINT_LINE);
        printf("[cco] Working on terminal number %d\n",this->num_terminals+1);
        fprintf(log_file,"%s\n",PRINT_LINE);
        fprintf(log_file,"[cco] Working on terminal number %d\n",this->num_terminals+1);

        sucess = generate_terminal(cost_function_config,local_opt_config); 

        // Write the current stage of the Purkinje network to a file
        if (sucess)
        {
            write_to_vtk_iteration();
        }

        printf("%s\n",PRINT_LINE);
        fprintf(log_file,"%s\n",PRINT_LINE);
    }
    // ===============================================================================================================

    // Check for the remaining active or inactive PMJ's to connect
    if (this->using_pmj_location)
    {
        printf("Number of connected PMJ's in the main loop: %u/%u\n",this->pmj_data->total_num_connected,this->pmj_data->connected.size());
        fprintf(log_file,"Number of connected PMJ's in the main loop: %u/%u\n",this->pmj_data->total_num_connected,this->pmj_data->connected.size());

        // Connect the remaining active PMJ's by loosening the LAT error threashold
        sucess = connect_remaining_active_pmjs(cost_function_config,local_opt_config);
        if (sucess) 
        {
            printf("[+] All active PMJ's were connected!\n");
        }
        else 
        {
            printf("[-] ERROR! Could not connect all PMJ's!\n");
            exit(EXIT_FAILURE);
        }

        // Connect the remaining inactives terminals
        sucess = connect_remaining_inactive_pmjs(cost_function_config,local_opt_config);

        if (sucess) write_to_vtk_iteration();
    }
}

void CCO_Network::make_root_using_cloud_points ()
{
    bool is_active;
    double x_prox[3], x_inew[3], ref_lat;
    
    uint32_t K_term = 1;
    double r_supp = R_PERF;

    // Set the position of the root
    memcpy(x_prox,this->root_pos,sizeof(double)*3);
    memset(x_inew,0.0,sizeof(double)*3);

    Point *A = new Point(0,x_prox);
    Point *B = new Point(1,x_inew);

    // Sort the distal position of the root using the cloud of points until the root segment has a size larger than a threashold distance
    uint32_t index = 0;
    uint32_t counter = 0;
    uint32_t iterations = 0;
    bool is_root_ok = false;
    double d_threashold = sqrt( M_PI * r_supp * r_supp );
    while (!is_root_ok)
    {
        uint32_t selected_index = this->cloud_data->sort_point(B,this->max_rand_offset);
        bool is_active = B->is_active;
        double root_length = euclidean_norm(A->x,A->y,A->z,B->x,B->y,B->z);

        //if (root_length >= d_threashold && !has_intersect_obstacle(x_prox,x_inew,obstacle_faces))
        if (root_length >= d_threashold && !is_active)
        {
            this->cloud_data->connected[selected_index] = true;
            is_root_ok = true;
        }
            
        else
            counter++;
        
        if (counter > 8) d_threashold *= 0.9;

        iterations++;   
    }
    printf("[cco] Root segment was set in %u iterations.\n", iterations);

    // Insert the root points in the point list
    this->point_list.push_back(A);
    this->point_list.push_back(B);

    // Create the root segment
    Segment *iroot = new Segment(0,A,B,NULL,NULL,NULL);
    this->segment_list.push_back(iroot);
    
    // Rescale the tree and update the number of terminals
    rescale_root(iroot);
    this->num_terminals = 1;

    // Write the current Purkinje tree in VTK
    write_to_vtk_iteration();
}

void CCO_Network::make_root_using_initial_network () { }

void CCO_Network::rescale_root (Segment *iroot)
{
    iroot->radius = this->start_radius;
}

void CCO_Network::rescale_tree (Segment *ibiff, Segment *iconn, Segment *inew)
{
    // [FRACTAL] Symmetric 
    double r_par = iconn->radius;
    double r_left = pow(0.5, 1.0/this->gamma) * r_par;
    double r_right = pow(0.5, 1.0/this->gamma) * r_par;

    // Fix a bifurcation ratio (Decrease the radius at each level of the tree by fixed factor)
    //inew->beta = r_left / r_par;
    //iconn->beta = r_right / r_par;

    // CONSTANT RADIUS
    double radius_ratio = r_left / r_right;
    inew->beta = radius_ratio;
    iconn->beta = radius_ratio;

    // Rescale the until we reach the root by using the "parent" pointer
    Segment *ipar = ibiff->parent;
    if (ipar != NULL)
    {
        Segment *ipar_left = ipar->left;
        Segment *ipar_right = ipar->right;
        rescale_until_root(ipar,ipar_left,ipar_right);
    }
    // We are already at the root
    else
    {
        ibiff->radius = this->start_radius;
    }
}

void CCO_Network::rescale_until_root (Segment *ipar, Segment *ipar_left, Segment *ipar_right)
{
    // Root was reached
    if (ipar == NULL) return;

    if (ipar_left != NULL && ipar_right != NULL)
    {
        // Symmetric
        double r_par = ipar->radius;
        double r_left = pow(0.5, 1.0/gamma) * r_par;
        double r_right = pow(0.5, 1.0/gamma) * r_par;
        //ipar_left->beta = r_left / r_par;
        //ipar_right->beta = r_right / r_par;
        
        // CONSTANT RADIUS
        double radius_ratio = r_left / r_right;
        ipar_left->beta = radius_ratio;
        ipar_right->beta = radius_ratio;

        // Call the function recursively until we reach the root
        if (ipar->parent != NULL)
            rescale_until_root(ipar->parent,ipar->parent->left,ipar->parent->right);
        else
            ipar->radius = this->start_radius;
    }
}

void CCO_Network::recalculate_radius ()
{
    for (uint32_t i = 0; i < this->segment_list.size(); i++)
    {
        Segment *cur_segment = this->segment_list[i];

        if (!cur_segment->adjusted)
            cur_segment->radius = cur_segment->calc_radius();
    }
}

void CCO_Network::recalculate_length ()
{
    for (uint32_t i = 0; i < this->segment_list.size(); i++)
    {
        Segment *cur_segment = this->segment_list[i];

        cur_segment->length = cur_segment->calc_length();
    }
}

void CCO_Network::restore_state_tree (Segment *iconn)
{
    Segment *ibiff = iconn->parent;
    Segment *inew = ibiff->left;
    Segment *ibiff_par = ibiff->parent;

    update_segment_pointers(iconn,ibiff,inew,true);

    // Eliminate "ibiff" and "inew"
    Point *M = inew->src;
    Point *T = inew->dest;
    eliminate_point_from_list(this->point_list,M);
    eliminate_point_from_list(this->point_list,T);
    eliminate_segment_from_list(this->segment_list,inew);
    eliminate_segment_from_list(this->segment_list,ibiff);

    // Update "ndist" from "iconn" until we reach the root
    this->num_terminals = update_ndist(iconn,true);

    // Update the 'betas' until we reach the root
    Segment *ipar = iconn->parent;
    if (ipar != NULL)
    {
        Segment *ipar_left = ipar->left;
        Segment *ipar_right = ipar->right;

        rescale_until_root(ipar,ipar_left,ipar_right);
    }
    else
    {
        iconn->radius = start_radius;
    }

    // Update the segments radius and length
    recalculate_radius();
    recalculate_length();
}

bool CCO_Network::generate_terminal_using_cloud_points (CostFunctionConfig *cost_function_config, LocalOptimizationConfig *local_opt_config)
{
    double new_pos[3];

    bool sucess = false;
    bool point_is_ok = false;
    uint32_t tosses = 0;
    uint32_t K_term = this->num_terminals;
    double d_threash = calc_dthreashold(R_PERF,K_term);     // (r_supp = R_PERF) in this case
    
    // Array of feasible segments to connect the new terminal
    std::vector<Segment*> feasible_segments;

    // Reference to the segment we are going to make the connection
    Segment *iconn = NULL;

    // Pointer to the new terminal point
    uint32_t cur_num_points = this->point_list.size();
    Point *new_term = new Point(cur_num_points+1);

    while (!point_is_ok)
    {
        // Reset the feasible segments list
        feasible_segments.clear();

        // Sort a terminal position from the cloud of points and update the 'new_term' data
        uint32_t selected_index = this->cloud_data->sort_point(new_term,this->max_rand_offset);

        // Check the distance criterion for this point
        point_is_ok = connection_search(this->segment_list,new_term,d_threash);

        if (point_is_ok) point_is_ok = check_collisions_and_fill_feasible_segments(this->segment_list,new_term,feasible_segments);

        if (point_is_ok)
        {
            // COST FUNCTION
            iconn = cost_fn->eval(this,cost_function_config,local_opt_config,\
                                feasible_segments,\
                                new_term,false);
            
            // No feasible point was found 
            if (iconn == NULL) point_is_ok = false;
        }

        // If the point does not attend the distance criterion or if there is a collision
        // we need to choose another point
        if (!point_is_ok)
        {
            tosses++;
            if (tosses > NTOSS)
            {
                fprintf(log_file,"[!] Reducing dthreash! Before = %g || Now = %g \n",\
                        d_threash,d_threash*FACTOR);
                d_threash *= FACTOR;
                tosses = 0;
            }

            // GIVE UP!
            if (d_threash < D_THREASH_LIMIT)
            {
                delete new_term;

                return sucess;
            }
                
        }
        // The point is a valid one and we can eliminate it from the cloud
        else
        {   
            printf("[!] Selected point %u -- (%g %g %g)\n",selected_index,\
                                                this->cloud_data->points[selected_index]->x,\
                                                this->cloud_data->points[selected_index]->y,\
                                                this->cloud_data->points[selected_index]->z);
            this->cloud_data->connected[selected_index] = true;
        }
    }

    Segment *inew = build_segment(local_opt_config,iconn->id,new_term);

    if (check_null_segments(this->segment_list))
    {
        sucess = false;
        restore_state_tree(iconn);
        exit(1);
    }
    else
    {
        sucess = true;
        update_min_max_terminal_lat();
        printf("Terminals --> min.LAT = %g || max.LAT = %g\n",this->min_term_lat,this->max_term_lat);
    }

    delete new_term;

    return sucess;
}

bool CCO_Network::generate_terminal_using_point (CostFunctionConfig *cost_function_config, LocalOptimizationConfig *local_opt_config, Point *pmj_point, const bool evaluate)
{
    bool sucess = false;
    bool point_is_ok = false;
    uint32_t counter_tries = 0;
    uint32_t K_term = this->num_terminals;
    double d_threash = calc_dthreashold(R_PERF,K_term);     
    
    // Array of feasible segments to connect the new terminal
    std::vector<Segment*> feasible_segments;

    // Reference to the segment we are going to make the connection
    Segment *iconn = NULL;

    while (!point_is_ok)
    {
        // Reset the feasible segments list
        feasible_segments.clear();

        // Check the distance criterion for this point
        point_is_ok = connection_search(this->segment_list,pmj_point,d_threash);

        if (point_is_ok) point_is_ok = check_collisions_and_fill_feasible_segments(this->segment_list,pmj_point,feasible_segments);

        if (point_is_ok)
        {
            // [COST FUNCTION]
            iconn = cost_fn->eval(this,cost_function_config,local_opt_config,\
                                feasible_segments,\
                                pmj_point,pmj_point->is_active);
            
            // No feasible point was found 
            if (iconn == NULL) point_is_ok = false;
        }

        // If the point does not attend the distance criterion or if there is a collision
        // we need to choose another point
        if (!point_is_ok)
        {
            fprintf(log_file,"[!] Reducing dthreash! Before = %g || Now = %g \n",d_threash,d_threash*FACTOR);
            d_threash *= FACTOR;
            counter_tries++;

            // GIVE UP!
            if (d_threash < D_THREASH_LIMIT)
            {
                return sucess;
            }    
        }
    }

    // The point is feasible, so we build the new segment of the tree
    Segment *inew = build_segment(local_opt_config,iconn->id,pmj_point);

    // [EVALUATE PMJ ACTIVATION]
    if (evaluate)
    {
        // Step 1: Connect using normal CCO
        sucess = evaluate_pmj_local_activation_time(inew,pmj_point,cost_function_config);
        if (sucess) printf("\t\t[+] PMJ point %u was connected in step 1!\n",pmj_point->id);

        // Step 2: Connect using normal CCO and region radius
        if (!sucess)
        {   
            sucess = attempt_connect_using_region_radius(pmj_point,cost_function_config);
            if (sucess) printf("\t\t[+] PMJ point %u was connected in step 2!\n",pmj_point->id);
        }

        point_is_ok = sucess;
        this->cloud_data->connected[pmj_point->id] = sucess;
        if (sucess) this->pmj_data->total_num_connected++;
    }
    else
    {
        sucess = true;
    }
    
    return sucess;
}

bool CCO_Network::generate_terminal (CostFunctionConfig *cost_function_config, LocalOptimizationConfig *local_opt_config)
{
    double new_pos[3];

    bool sucess = false;
    bool point_is_ok = false;
    uint32_t counter_tries = 0;
    uint32_t tosses = 0;
    uint32_t K_term = this->num_terminals;
    double d_threash = calc_dthreashold(R_PERF,K_term);     // (r_supp = R_PERF) in this case
    
    // Array of feasible segments to connect the new terminal
    std::vector<Segment*> feasible_segments;

    // Reference to the segment we are going to make the connection
    Segment *iconn = NULL;

    // Pointer to the new terminal point
    uint32_t cur_num_points = this->point_list.size();
    Point *new_term = new Point(cur_num_points+1);

    while (!point_is_ok)
    {
        // Reset the feasible segments list
        feasible_segments.clear();

        // Sort a terminal position from the cloud of points and update the 'new_term' data
        uint32_t selected_index = this->cloud_data->sort_point(new_term,this->max_rand_offset);

        // Check the distance criterion for this point
        point_is_ok = connection_search(this->segment_list,new_term,d_threash);

        if (point_is_ok) point_is_ok = check_collisions_and_fill_feasible_segments(this->segment_list,new_term,feasible_segments);

        if (point_is_ok)
        {
            // COST FUNCTION
            iconn = cost_fn->eval(this,cost_function_config,local_opt_config,\
                                feasible_segments,\
                                new_term,new_term->is_active);
            
            // No feasible point was found 
            if (iconn == NULL) point_is_ok = false;
        }

        // If the point does not attend the distance criterion or if there is a collision
        // we need to choose another point
        if (!point_is_ok)
        {
            tosses++;
            if (tosses > NTOSS)
            {
                fprintf(log_file,"[!] Reducing dthreash! Before = %g || Now = %g \n",\
                        d_threash,d_threash*FACTOR);
                d_threash *= FACTOR;
                counter_tries++;
                tosses = 0;
            }

            // GIVE UP!
            if (d_threash < D_THREASH_LIMIT)
            {
                delete new_term;

                return sucess;
            }
                
        }
        // The point is a valid one and we can eliminate it from the cloud
        else
        {   
            printf("[!] Selected point %u -- (%g %g %g) {Active=%d}\n",selected_index,\
                                                this->cloud_data->points[selected_index]->x,\
                                                this->cloud_data->points[selected_index]->y,\
                                                this->cloud_data->points[selected_index]->z,\
                                                this->cloud_data->points[selected_index]->is_active);

            Segment *inew = build_segment(local_opt_config,iconn->id,new_term);

            // [EVALUATE PMJ ACTIVATION]
            if (new_term->is_active)
            {
                // Step 1: Connect using normal CCO
                sucess = evaluate_pmj_local_activation_time(inew,new_term,cost_function_config);
                if (sucess) printf("\t\t[+] PMJ point %u was connected in step 1!\n",new_term->id);

                // Step 2: Connect using normal CCO and region radius
                if (!sucess)
                {   
                    sucess = attempt_connect_using_region_radius(new_term,cost_function_config);
                    if (sucess) printf("\t\t[+] PMJ point %u was connected in step 2!\n",new_term->id);
                }
                
                point_is_ok = sucess;
                this->cloud_data->connected[selected_index] = sucess;
                if (sucess) this->pmj_data->total_num_connected++;
            }
            // It is not an active PMJ, so just connect and leave
            else
            {
                sucess = true;
            }
        }
    }

    update_min_max_terminal_lat();
    printf("Terminals --> min.LAT = %g || max.LAT = %g\n",this->min_term_lat,this->max_term_lat);

    delete new_term;

    return sucess;
}

Segment* CCO_Network::build_segment (LocalOptimizationConfig *local_opt_config, const uint32_t index, Point *new_term)
{
    uint32_t cur_num_segments = this->segment_list.size();
    Segment *iconn = this->segment_list[index];
    
    // Create the bifurcation point
    Point *M = generate_bifurcation_node(this->point_list,iconn,local_opt_config,this->using_local_optimization);
    this->point_list.push_back(M);

    // Create the terminal point
    Point *T = generate_terminal_node(this->point_list,new_term);
    this->point_list.push_back(T);

    // Create 'ibiff'
    Segment *ibiff = new Segment(cur_num_segments,iconn->src,M,\
                                NULL,NULL,iconn->parent);
    ibiff->ndist = iconn->ndist;
    this->segment_list.push_back(ibiff);

    // Create 'inew'
    Segment *inew = new Segment(cur_num_segments+1,M,T,\
                                NULL,NULL,ibiff);
    this->segment_list.push_back(inew);

    // Update pointers
    update_segment_pointers(iconn,ibiff,inew,false);

    // Update ndist from ibiff until the root
    this->num_terminals = update_ndist(ibiff,false);

    // Rescale the network from the segments of the new bifurcation
    rescale_tree(ibiff,iconn,inew);

    // Recalculate the segment radius using the 'beta' values
    recalculate_radius();
    recalculate_length();

    // Return a reference to the new node
    return inew;
}

void CCO_Network::prune_segment (Segment *inew)
{
    Segment *iconn = NULL;
    Segment *ibiff_par = NULL;
    Segment *ibiff = NULL;

    // Set the pointers 
    ibiff = inew->parent;
    if (ibiff->left == inew) iconn = ibiff->right;
    else if (ibiff->right == inew) iconn = ibiff->left;
    ibiff_par = ibiff->parent;

    // Update "iconn" pointers and values
    iconn->parent = ibiff->parent;
    iconn->src = ibiff->src;
    iconn->beta = ibiff->beta;

    // Update "ibiff" parent subtree pointer if exists
    if (ibiff_par != NULL)
    {
        if (ibiff_par->right == ibiff) ibiff_par->right = iconn;
        else if (ibiff_par->left == ibiff) ibiff_par->left = iconn;
    }

    // Eliminate the bifurcation and terminal point
    Point *M = inew->src;
    Point *T = inew->dest;
    eliminate_point_from_list(this->point_list,M);
    eliminate_point_from_list(this->point_list,T);
    eliminate_segment_from_list(this->segment_list,inew);
    eliminate_segment_from_list(this->segment_list,ibiff);

    // Update "ndist" from "iconn" until we reach the root
    this->num_terminals = update_ndist(iconn,true);

    // Update the 'betas' until we reach the root
    Segment *ipar = iconn->parent;
    if (ipar != NULL)
    {
        Segment *ipar_left = ipar->left;
        Segment *ipar_right = ipar->right;

        rescale_until_root(ipar,ipar_left,ipar_right);
    }
    else
    {
        iconn->radius = start_radius;
    }

    // Update the segments radius and length
    recalculate_radius();
    recalculate_length();
}

bool CCO_Network::evaluate_pmj_local_activation_time (Segment *inew, Point *pmj_point, CostFunctionConfig *cost_function_config)
{
    double lat = inew->calc_terminal_local_activation_time() + this->lat_offset;
    double lat_error = (pmj_point->lat - lat);

    // Check if the LAT error is less then the absolute error tolerance
    if (fabs(lat_error) < this->pmj_data->lat_error_tolerance)
    {
        // Try to adjust the diameter to improve even more the error
        bool sucess = adjust_terminal_diameter(inew,pmj_point);
        if (sucess)
        {
            lat = inew->calc_terminal_local_activation_time() + this->lat_offset;
            lat_error = (pmj_point->lat - lat);
        }

        inew->can_touch = false;
        this->cloud_data->error[pmj_point->id] = lat_error;
        this->cloud_data->aprox[pmj_point->id] = lat;
        printf("\t[PMJ %u] Ref: %g ms || Aprox: %g ms || Error: %g ms\n",pmj_point->id,pmj_point->lat,lat,lat_error);
        return true;
    }
    // If the error is above the tolerance we will try adjust the diameter
    else
    {
        bool sucess = adjust_terminal_diameter(inew,pmj_point);
        if (sucess)
        {
            lat = inew->calc_terminal_local_activation_time() + this->lat_offset;
            lat_error = (pmj_point->lat - lat);
            inew->can_touch = false;
            this->cloud_data->error[pmj_point->id] = lat_error;
            this->cloud_data->aprox[pmj_point->id] = lat;
            printf("\t[PMJ %u] Ref: %g ms || Aprox: %g ms || Error: %g ms\n",pmj_point->id,pmj_point->lat,lat,lat_error);
        }
        else
        {
            prune_segment(inew);
        }
        return false;
    }
}

bool CCO_Network::adjust_terminal_diameter (Segment *inew, Point *pmj_point)
{
    double original_radius = inew->radius;
    double original_velocity = inew->calc_propagation_velocity();

    double lat = inew->calc_terminal_local_activation_time() + this->lat_offset;
    double lat_error = (pmj_point->lat - lat);

    // The propagation velocity is too fast for that PMJ
    if (lat_error > 0.0)
    {
        // Try to adjust the diameter of the last segment by decreasing its velocity
        double lat_parent = inew->parent->calc_terminal_local_activation_time() + this->lat_offset;
        double delta_t_target = pmj_point->lat - lat_parent;
        
        // If the PMJ has been activated too soon we adjust the diameter in order to fit the reference LAT
        if (delta_t_target > 0.0)
        {
            double delta_s_term = inew->calc_pathway_length()*M_TO_UM;
            double delta_s_parent = inew->parent->calc_pathway_length()*M_TO_UM;
            double new_cv = (delta_s_term - delta_s_parent) / delta_t_target;
            
            // We only update the segment diameter when the new CV is greater than 1000um/ms --> {1m/s}
            if (new_cv > MIN_CV_THREASHOLD && new_cv < MAX_CV_THREASHOLD)
            {
                double prev_error = lat_error;
                inew->update_radius(new_cv);
            
                lat = inew->calc_terminal_local_activation_time() + this->lat_offset;
                lat_error = (pmj_point->lat - lat);

                printf("\t[PMJ %u] Adjusting diameter from %g um to %g um || Propagation velocity decreased from %g to %g m/s || Error decreased from %g to %g ms\n",pmj_point->id,original_radius*2.0*MM_TO_UM,inew->radius*2.0*MM_TO_UM,original_velocity,new_cv/1000,prev_error,lat_error);
                printf("\t[PMJ %u] Ref: %g ms || Aprox: %g ms || Error: %g ms\n",pmj_point->id,pmj_point->lat,lat,lat_error);

                if (fabs(lat_error) < this->pmj_data->lat_error_tolerance)
                {
                    inew->can_touch = false;
                    this->cloud_data->error[pmj_point->id] = lat_error;
                    this->cloud_data->aprox[pmj_point->id] = lat;
                    return true;
                }
            }
        }
    }
    return false;
}

bool CCO_Network::attempt_connect_using_region_radius (Point *pmj_point, CostFunctionConfig *cost_function_config)
{
    bool sucess = false;
    double center[3], ori_pos[3], original_radius, original_velocity, lat_error, lat;
    double radius = this->pmj_data->region_radius;

    center[0] = pmj_point->x;
    center[1] = pmj_point->y;
    center[2] = pmj_point->z;

    for (uint32_t i = 0; i < this->segment_list.size(); i++)
    {
        Segment *cur_segment = this->segment_list[i];
        if (cur_segment->is_terminal() && cur_segment->is_inside_region(center,radius) && cur_segment->is_touchable())
        {
            // Save the original position and radius
            ori_pos[0] = cur_segment->dest->x;
            ori_pos[1] = cur_segment->dest->y;
            ori_pos[2] = cur_segment->dest->z;
            original_radius = cur_segment->radius;
            original_velocity = cur_segment->calc_propagation_velocity();

            // Change coordinate position
            cur_segment->dest->x = pmj_point->x;
            cur_segment->dest->y = pmj_point->y;
            cur_segment->dest->z = pmj_point->z;
            cur_segment->length = cur_segment->calc_length();
            cur_segment->can_touch = false;

            lat = cur_segment->calc_terminal_local_activation_time() + this->lat_offset;
            lat_error = (lat - pmj_point->lat);

            // Check if the LAT error is less then the absolute error tolerance
            if (fabs(lat_error) < this->pmj_data->lat_error_tolerance)
            {
                // Try to adjust the diameter to improve even more the error
                bool sucess = adjust_terminal_diameter(cur_segment,pmj_point);
                if (sucess)
                {
                    lat = cur_segment->calc_terminal_local_activation_time() + this->lat_offset;
                    lat_error = (pmj_point->lat - lat);
                }

                printf("\t[PMJ %u] Ref: %g ms || Aprox: %g ms || Error: %g ms\n",pmj_point->id,pmj_point->lat,lat,lat_error);
                cur_segment->can_touch = false;
                this->cloud_data->error[pmj_point->id] = lat_error;
                this->cloud_data->aprox[pmj_point->id] = lat;
                return true;
            }
            else
            {
                cur_segment->dest->x = ori_pos[0];
                cur_segment->dest->y = ori_pos[1];
                cur_segment->dest->z = ori_pos[2];
                cur_segment->length = cur_segment->calc_length();
                cur_segment->radius = original_radius;
                cur_segment->can_touch = true;
            }
        }
    }

    return false;
}

void CCO_Network::write_to_vtk_iteration ()
{
    uint32_t num_points = this->point_list.size();
    uint32_t num_segments = this->segment_list.size();
    uint32_t num_terminals = this->num_terminals;
    
    std::string filename;
    std::ostringstream os;
    os << this->output_dir << "/tree_nterm_" << num_terminals << ".vtk";
    filename = os.str();

    FILE *file = fopen(filename.c_str(),"w+");

    // Write the header
    fprintf(file,"# vtk DataFile Version 3.0\n");
    fprintf(file,"Tree\n");
    fprintf(file,"ASCII\n");
    fprintf(file,"DATASET POLYDATA\n");
    fprintf(file,"POINTS %u float\n",num_points);
    for (uint32_t i = 0; i < num_points; i++)
    {
        fprintf(file,"%g %g %g\n",this->point_list[i]->x * M_TO_UM,\
                                this->point_list[i]->y * M_TO_UM,\
                                this->point_list[i]->z * M_TO_UM);
    }
    
    fprintf(file,"LINES %u %u\n",num_segments,num_segments*3);
    for (uint32_t i = 0; i < num_segments; i++)
    {
        fprintf(file,"2 %u %u\n",this->segment_list[i]->src->id,\
                                this->segment_list[i]->dest->id);
    }
    
    fprintf(file,"CELL_DATA %u\n",num_segments);
    fprintf(file,"SCALARS radius float\n");
    fprintf(file,"LOOKUP_TABLE default\n");
    for (uint32_t i = 0; i < num_segments; i++)
    {
        fprintf(file,"%g\n",this->segment_list[i]->radius * MM_TO_UM);
    }

    fclose(file);
}

CCO_Network* CCO_Network::copy ()
{
    CCO_Network *result = new CCO_Network();

    result->num_terminals = this->num_terminals;
    result->N_term = this->N_term;
    result->seed = this->seed;
    result->max_rand_offset = this->max_rand_offset;
    result->gamma = this->gamma;
    result->lat_offset = this->lat_offset;
    result->cur_rand_index = 0;
    result->cur_pmj_package = 0;
    result->max_lat_error = __DBL_MIN__;
    result->min_max_aprox_lat[0] = __DBL_MAX__; result->min_max_aprox_lat[1] = __DBL_MIN__;
    result->min_max_ref_lat[0] = __DBL_MAX__; result->min_max_ref_lat[1] = __DBL_MIN__;

    memcpy(result->root_pos,this->root_pos,sizeof(double)*3);

    result->using_only_murray_law = this->using_only_murray_law;
    result->start_radius = this->start_radius;
    
    result->using_cloud_points = this->using_cloud_points;
    result->using_local_optimization = this->using_local_optimization;
    result->using_pmj_location = this->using_pmj_location;
        
    memcpy(result->min_max_aprox_lat,this->min_max_aprox_lat,sizeof(double)*2);
    memcpy(result->min_max_ref_lat,this->min_max_ref_lat,sizeof(double)*2);
    
    result->output_dir = "outputs/copy";
    result->cost_function_name = this->cost_function_name;
    result->local_optimization_function_name = this->local_optimization_function_name;
    
    result->Q_TERM = Q_PERF / result->N_term;
    result->R_PERF = powf(((3.0*V_PERF)/(4.0 * M_PI)), (1.0/3.0));

    result->set_cost_function();
    result->set_local_optimization_function_name();

    // Copy the Point array
    for (uint32_t i = 0; i < this->point_list.size(); i++)
    {
        Point *p = new Point(this->point_list[i]);
        result->point_list.push_back(p);
    }

    // Copy the Segment array
    for (uint32_t i = 0; i < this->segment_list.size(); i++)
    {
        uint32_t src_index = this->segment_list[i]->src->id;
        uint32_t dest_index = this->segment_list[i]->dest->id;
        Point *src = result->point_list[src_index];
        Point *dest = result->point_list[dest_index];

        Segment *s = new Segment(i,src,dest,NULL,NULL,NULL);
        s->ndist = this->segment_list[i]->ndist;
        s->beta = this->segment_list[i]->beta;
        s->radius = this->segment_list[i]->radius;
        s->prune = this->segment_list[i]->prune;
        s->can_touch = this->segment_list[i]->can_touch;

        result->segment_list.push_back(s);
    }
    // Set the Segment pointers
    for (uint32_t i = 0; i < this->segment_list.size(); i++)
    {
        Segment *left = this->segment_list[i]->left;
        Segment *right = this->segment_list[i]->right;
        Segment *parent = this->segment_list[i]->parent;

        result->segment_list[i]->left = (left) ? result->segment_list[left->id] : NULL;
        result->segment_list[i]->right = (right) ? result->segment_list[right->id] : NULL;
        result->segment_list[i]->parent = (parent) ? result->segment_list[parent->id] : NULL;
    }

    // Copy the Cloud data
    if (this->cloud_data)
    {
        result->cloud_data = this->cloud_data->copy();
    }
    
    // Copy the PMJ data
    if (this->pmj_data)
    {
        result->pmj_data = this->pmj_data->copy();
    }

    return result;
}

CCO_Network* CCO_Network::concatenate (CCO_Network *input)
{
    uint32_t offset_points, offset_segments;

    // Get a copy of the current network
    CCO_Network *result = this->copy();

    // Sum the number of terminals from both networks and update the PMJ flag
    result->num_terminals += input->num_terminals;
    result->using_pmj_location |= input->using_pmj_location;

    // Adjust the Min/Max PMJ errors
    result->min_max_aprox_lat[0] = (input->min_max_aprox_lat[0] < this->min_max_aprox_lat[0]) ? input->min_max_aprox_lat[0] : this->min_max_aprox_lat[0];
    result->min_max_aprox_lat[1] = (input->min_max_aprox_lat[1] > this->min_max_aprox_lat[1]) ? input->min_max_aprox_lat[1] : this->min_max_aprox_lat[1];
    result->min_max_ref_lat[0] = (input->min_max_ref_lat[0] < this->min_max_ref_lat[0]) ? input->min_max_ref_lat[0] : this->min_max_ref_lat[0];
    result->min_max_ref_lat[1] = (input->min_max_ref_lat[1] > this->min_max_ref_lat[1]) ? input->min_max_ref_lat[1] : this->min_max_ref_lat[1];
    result->max_lat_error = (input->max_lat_error > this->max_lat_error) ? input->max_lat_error : this->max_lat_error;

    // Concatenate the points
    offset_points = result->point_list.size();
    for (uint32_t i = 0; i < input->point_list.size(); i++)
    {
        Point *p = new Point(input->point_list[i]);
        p->id += offset_points;
        result->point_list.push_back(p);
    }

    // Concatenate the segments
    offset_segments = result->segment_list.size();
    for (uint32_t i = 0; i < input->segment_list.size(); i++)
    {
        uint32_t src_index = input->segment_list[i]->src->id + offset_points;
        uint32_t dest_index = input->segment_list[i]->dest->id + offset_points;
        Point *src = result->point_list[src_index];
        Point *dest = result->point_list[dest_index];

        Segment *s = new Segment(i,src,dest,NULL,NULL,NULL);
        s->id += offset_segments;
        s->ndist = input->segment_list[i]->ndist;
        s->beta = input->segment_list[i]->beta;
        s->radius = input->segment_list[i]->radius;
        s->prune = input->segment_list[i]->prune;
        s->can_touch = input->segment_list[i]->can_touch;

        result->segment_list.push_back(s);
    }
    // Update the segment pointers
    for (uint32_t i = 0; i < input->segment_list.size(); i++)
    {
        Segment *left = input->segment_list[i]->left;
        Segment *right = input->segment_list[i]->right;
        Segment *parent = input->segment_list[i]->parent;

        result->segment_list[i + offset_segments]->left = (left) ? result->segment_list[left->id + offset_segments] : NULL;
        result->segment_list[i + offset_segments]->right = (right) ? result->segment_list[right->id + offset_segments] : NULL;
        result->segment_list[i + offset_segments]->parent = (parent) ? result->segment_list[parent->id + offset_segments] : NULL;
    }

    // Concatenate cloud of points
    result->cloud_data->concatenate(input->cloud_data);

    // Concatenate PMJ's
    if (input->pmj_data)
    {
        result->pmj_data->concatenate(input->pmj_data);
    }

    return result;
}

void CCO_Network::link_segments (Segment *term_1, Segment *term_2)
{
    // New version
    double middle_pos[3];
    uint32_t cur_num_points = this->point_list.size();
    uint32_t cur_num_segments = this->segment_list.size();

    term_1->calc_middle_point(middle_pos);
    Point *M = new Point(cur_num_points,middle_pos);
    this->point_list.push_back(M);

    Segment *ibiff = term_1;
    Segment *iconn = new Segment(cur_num_segments,M,term_1->dest,NULL,NULL,NULL);
    Segment *inew = new Segment(cur_num_segments+1,M,term_2->src,NULL,NULL,NULL);

    // ibiff
    ibiff->dest = M;
    ibiff->length = ibiff->calc_length();
    ibiff->right = inew;
    ibiff->left = iconn;

    // iconn
    iconn->parent = ibiff;

    // inew 
    inew->parent = ibiff;
    inew->right = term_2;

    this->segment_list.push_back(iconn);
    this->segment_list.push_back(inew);
}

void CCO_Network::adjust_radius ()
{
    // Set the output directory first
    set_save_network();

    for (uint32_t i = 0; i < this->segment_list.size(); i++)
    {
        Segment *cur_segment = this->segment_list[i];

        uint32_t level = cur_segment->calc_level()-1;

        double radius = this->start_radius;
        for (uint32_t j = 0; j < level; j++)
            radius *= pow(0.5, 1.0/this->gamma);
        
        cur_segment->radius = radius;
    }    

    // TODO: Think about we are going adjust the 'ndist' and 'beta' variables
}

void CCO_Network::adjust_radius_2 ()
{
    // Set the output directory first
    set_save_network();

    for (uint32_t i = 0; i < this->segment_list.size(); i++)
    {
        Segment *cur_segment = this->segment_list[i];

        cur_segment->radius = this->start_radius;
    }    
}

void CCO_Network::calc_electric_error ()
{
    // Compute min/max LAT values from the active PMJ's (aproximation)
    for (uint32_t i = 0; i < this->cloud_data->aprox.size(); i++)
    {
        uint32_t id = this->cloud_data->points[i]->id;
        if (this->cloud_data->connected[id] && this->cloud_data->points[id]->is_active)
        {
            double lat = this->cloud_data->aprox[id];
            if (lat < this->min_max_aprox_lat[0]) this->min_max_aprox_lat[0] = lat;
            if (lat > this->min_max_aprox_lat[1]) this->min_max_aprox_lat[1] = lat;
        }  
    }

    // Compute min/max LAT values from the active PMJ's (reference)
    for (uint32_t i = 0; i < this->cloud_data->points.size(); i++)
    {
        uint32_t id = this->cloud_data->points[i]->id;
        if (this->cloud_data->connected[id] && this->cloud_data->points[id]->is_active)
        {
            double lat = this->cloud_data->points[id]->lat;
            if (lat < this->min_max_ref_lat[0]) this->min_max_ref_lat[0] = lat;
            if (lat > this->min_max_ref_lat[1]) this->min_max_ref_lat[1] = lat;
        }
    }

    // Compute max LAT error from the active PMJ's
    for (uint32_t i = 0; i < this->cloud_data->error.size(); i++)
    {
        uint32_t id = this->cloud_data->points[i]->id;
        if (this->cloud_data->connected[id] && this->cloud_data->points[id]->is_active)
        {
            double error = fabs(this->cloud_data->error[id]);
            if (error > this->max_lat_error) this->max_lat_error = error;
        }
    }

    // Compute the RMSE and RRMSE from the active PMJ's
    uint32_t n = this->cloud_data->points.size();
    uint32_t counter = 0;
    double sum_num = 0.0;
    double sum_den = 0.0;
    for (uint32_t i = 0; i < n; i++)
    {
        uint32_t id = this->cloud_data->points[i]->id;
        if (this->cloud_data->connected[id] && this->cloud_data->points[id]->is_active)
        {
            double ref_value = this->cloud_data->points[id]->lat;
            double error = fabs(this->cloud_data->error[id]);

            sum_num += powf(error,2);
            sum_den += powf(ref_value,2);
            counter++;
        }
    }    
    double l2_norm = sqrt(sum_den);
    this->rmse = sqrt(sum_num/(double)counter);
    this->rrmse = sqrt(sum_num/sum_den);

    // Compute the number of active PMJ's that have an error less than a certain threashold
    uint32_t counter_less_2ms = 0;
    uint32_t counter_less_5ms = 0;
    for (uint32_t i = 0; i < this->cloud_data->error.size(); i++)
    {
        uint32_t id = this->cloud_data->points[i]->id;
        if (this->cloud_data->connected[id] && this->cloud_data->points[id]->is_active)
        {
            double error = fabs(this->cloud_data->error[id]);

            if (error < 2.0)    counter_less_2ms++;
            if (error < 5.0)    counter_less_5ms++;
        }
        
    }
    this->epsilon_2ms = (double)counter_less_2ms / (double)counter;
    this->epsilon_5ms = (double)counter_less_5ms / (double)counter;
}


bool CCO_Network::force_connection_to_closest_segment (CostFunctionConfig *cost_function_config, LocalOptimizationConfig *local_opt_config, Point *pmj)
{
    Segment *closest = NULL;
    double min_dist = __DBL_MAX__;

    for (uint32_t i = 0; i < this->segment_list.size(); i++)
    {
        Segment *cur_segment = this->segment_list[i];

        if (cur_segment->is_terminal() && cur_segment->is_touchable())
        {
            double dist = euclidean_norm(cur_segment->dest->x,cur_segment->dest->y,cur_segment->dest->z,\
                                        pmj->x,pmj->y,pmj->z);
            if (dist < min_dist)
            {
                min_dist = dist;
                closest = cur_segment;
            }
        }
    }

    if (closest)
    {
        // Change the coordinates
        closest->dest->x = pmj->x;
        closest->dest->y = pmj->y;
        closest->dest->z = pmj->z;
        closest->length = closest->calc_length();
        closest->can_touch = false;
        
        // Calculate the LAT error
        double lat = closest->calc_terminal_local_activation_time() + this->lat_offset;
        double lat_error = fabs(lat-pmj->lat);
        this->cloud_data->error[pmj->id] = lat_error;
        this->cloud_data->aprox[pmj->id] = lat;
        this->cloud_data->connected[pmj->id] = true;
        this->pmj_data->total_num_connected++;
        printf("\t[PMJ %u] Ref: %g ms || Aprox: %g ms || Error: %g ms\n",pmj->id,pmj->lat,lat,lat_error);
        printf("\t\tPMJ point %u was connected by force!\n",pmj->id);

        return true;
    }
    else
    {
        // Try to connect with the closest segment that is not a terminal by making a bifurcation ...

        printf("[-] ERROR! Could not force PMJ connection with PMJ point %u!\n",pmj->id);
        exit(EXIT_FAILURE);
    }
    return false;
}

bool CCO_Network::connect_remaining_active_pmjs (CostFunctionConfig *cost_function_config, LocalOptimizationConfig *local_opt_config)
{    
    bool sucess;
    double original_lat_tolerance = this->pmj_data->lat_error_tolerance;
    for (uint32_t i = 0; i < this->cloud_data->points.size(); i++)
    {
        bool is_connected = this->cloud_data->connected[i];
        bool is_active = this->cloud_data->points[i]->is_active;
        Point *pmj_point = this->cloud_data->points[i];
        sucess = false;

        // Check if the point is a non-connected active PMJ
        if (!is_connected && is_active)
        {
            // Until the PMJ point is not connected we will not leave this loop
            while (!sucess)
            {
                // Loose the LAT error tolerance by the PMJ_LOOSE_RATE
                printf("[!] Loosing PMJ error tolerance %g\n",this->pmj_data->lat_error_tolerance);
                this->pmj_data->lat_error_tolerance *= PMJ_LOOSE_RATE;

                sucess = generate_terminal_using_point(cost_function_config,local_opt_config,pmj_point,true);

                // If we surpass a certain loose limit we force connection to the closest segment in tree
                if (!sucess && this->pmj_data->lat_error_tolerance > PMJ_LOOSE_LIMIT)
                {
                    sucess = force_connection_to_closest_segment(cost_function_config,local_opt_config,pmj_point);
                }
            }
            // Reset the LAT error tolerance
            this->pmj_data->lat_error_tolerance = original_lat_tolerance;

        }

        if (sucess) write_to_vtk_iteration();
    }

    // Check if all active PMJ's are connected
    return check_active_pmjs_connection();
}

bool CCO_Network::connect_remaining_inactive_pmjs (CostFunctionConfig *cost_function_config, LocalOptimizationConfig *local_opt_config)
{
    bool sucess;
    uint32_t num_connected_pmjs = this->pmj_data->total_num_connected;
    uint32_t num_remanaing_inactives = (this->N_term + num_connected_pmjs) - (this->num_terminals);
    for (uint32_t i = 0; i < num_remanaing_inactives; i++)
    {
        printf("%s\n",PRINT_LINE);
        printf("[cco] Working on inactive terminal number %d\n",this->num_terminals+1);
        fprintf(log_file,"%s\n",PRINT_LINE);
        fprintf(log_file,"[cco] Working on inactive terminal number %d\n",this->num_terminals+1);

        sucess = generate_terminal(cost_function_config,local_opt_config);

        printf("%s\n",PRINT_LINE);
        fprintf(log_file,"%s\n",PRINT_LINE);
    }
    return sucess;
}

bool CCO_Network::attempt_connect_using_inverse (CostFunctionConfig *cost_function_config, LocalOptimizationConfig *local_opt_config, std::vector<Segment*> feasible_segments, Point *pmj_point)
{
    bool sucess = false;
    
    // [COST FUNCTION] = Maximize the total length
    Segment *iconn = this->pmj_data->cost_fn->eval(this,cost_function_config,local_opt_config,\
                        feasible_segments,\
                        pmj_point,true);
            
    // No feasible point was found 
    if (iconn == NULL) return sucess;

    // Build the new segment of the tree
    Segment *inew = build_segment(local_opt_config,iconn->id,pmj_point);

    // [EVALUATE PMJ ACTIVATION]
    // Step 3: Connect using inverse CCO
    sucess = evaluate_pmj_local_activation_time(inew,pmj_point,cost_function_config);
    
    return sucess;    
}

void CCO_Network::print ()
{
    print_point_list();
    print_segment_list();

    printf("Q_PERF = %g\n",Q_PERF);
    printf("Q_TERM = %g\n",Q_TERM);
    printf("P_PERF = %g\n",P_PERF);
    printf("P_TERM = %g\n",P_TERM);
    printf("V_PERF = %g\n",V_PERF);
    printf("R_PERF = %g\n",R_PERF);
}

void CCO_Network::print_point_list ()
{
    for (uint32_t i = 0; i < this->point_list.size(); i++)
        this->point_list[i]->print();
    printf("Total number of points = %u\n",this->point_list.size());
}

void CCO_Network::print_segment_list ()
{
    for (uint32_t i = 0; i < this->segment_list.size(); i++)
        this->segment_list[i]->print();
    printf("Total number of segments = %u\n",this->segment_list.size());
}

void CCO_Network::print_network_info ()
{
    std::vector<double> segments;
    get_segment_length(this->segment_list,segments);

    std::vector<double> angles;
    get_bifurcation_angles(this->segment_list,angles);

    double mean_segment_length, std_segment_length;
    calc_mean_std(segments,mean_segment_length,std_segment_length);
    write_vector_to_file(segments,this->output_dir + "/segments_length.dat");

    double mean_biff_angle, std_biff_angle;
    calc_mean_std(angles,mean_biff_angle,std_biff_angle);
    write_vector_to_file(angles,this->output_dir + "/bifurcation_angle.dat");

    write_geometric_info_to_file(this->output_dir + "/network_info.txt",\
                    segments,mean_segment_length,std_segment_length,\
                    angles,mean_biff_angle,std_biff_angle);
    //write_geometric_info_for_excel(this->output_dir + "/network_info_for_excel.txt",\
                    segments,mean_segment_length,std_segment_length,\
                    angles,mean_biff_angle,std_biff_angle);
    //write_geometric_info_for_latex(this->output_dir + "/network_info_for_latex.txt",\
                    segments,mean_segment_length,std_segment_length,\
                    angles,mean_biff_angle,std_biff_angle);
    
    printf("[INFO] Total number of segment = %u\n",segments.size());
    printf("[INFO] Segment length = %.2lf +/- %.2lf mm\n",mean_segment_length,std_segment_length);
    printf("[INFO] Total number of bifurcations = %u\n",angles.size());
    printf("[INFO] Bifurcation angle = %.2lf +/- %.2lf degrees\n",mean_biff_angle,std_biff_angle);
    if (this->using_pmj_location) 
    {
        printf("[INFO] Number of PMJ's connected = %u/%u\n",this->pmj_data->total_num_connected,this->pmj_data->points.size());
        write_vector_to_file(this->cloud_data->error,this->output_dir + "/cloud_error.dat");

        calc_electric_error();
        printf("[INFO] Max Error = %.2lf ms || Ref Min LAT = %.2lf ms || Ref Max LAT = %.2lf ms ||\n",this->max_lat_error,this->min_max_ref_lat[0],this->min_max_ref_lat[1]);                                                                                                           
        printf("                         || Aprox Min LAT = %.2lf ms || Aprox Max LAT = %.2lf ms ||\n",this->min_max_aprox_lat[0],this->min_max_aprox_lat[1]);
        printf("       RMSE = %.2lf ms || RRMSE = %.2lf %%\n",this->rmse,this->rrmse*100.0);    
        printf("       Epsilon 2ms = %.2lf %% || Epsilon = 5ms = %.2lf %%\n",this->epsilon_2ms*100.0,this->epsilon_5ms*100.0);       

        write_electric_info_to_file(this->output_dir + "/network_info.txt",\
                    this->max_lat_error,this->min_max_ref_lat[0],this->min_max_ref_lat[1],\
                    this->min_max_aprox_lat[0],this->min_max_aprox_lat[1],\
                    this->rmse,this->rrmse*100.0,\
                    this->epsilon_2ms*100.0,this->epsilon_5ms*100.0);
        //write_electric_info_for_excel(this->output_dir + "/network_info_for_excel.txt",\
                    this->max_lat_error,this->min_max_ref_lat[0],this->min_max_ref_lat[1],\
                    this->min_max_aprox_lat[0],this->min_max_aprox_lat[1],\
                    this->rmse,this->rrmse*100.0,\
                    this->epsilon_2ms*100.0,this->epsilon_5ms*100.0);
        //write_electric_info_for_latex(this->output_dir + "/network_info_for_latex.txt",\
                    this->max_lat_error,this->min_max_ref_lat[0],this->min_max_ref_lat[1],\
                    this->min_max_aprox_lat[0],this->min_max_aprox_lat[1],\
                    this->rmse,this->rrmse*100.0,\
                    this->epsilon_2ms*100.0,this->epsilon_5ms*100.0);                                                                                                 
    }
}

void CCO_Network::update_min_max_terminal_lat ()
{
    for (uint32_t i = 0; i < this->segment_list.size(); i++)
    {
        Segment *cur_segment = this->segment_list[i];
        if (cur_segment->is_terminal())
        {
            double lat = cur_segment->calc_terminal_local_activation_time();
            if (lat < min_term_lat) min_term_lat = lat;
            if (lat > max_term_lat) max_term_lat = lat;
        }
    }
}

bool CCO_Network::check_active_pmjs_connection ()
{
    for (uint32_t i = 0; i < this->cloud_data->points.size(); i++)
    {
        bool is_connected = this->cloud_data->connected[i];
        bool is_active = this->cloud_data->points[i]->is_active;

        if (is_active && !is_connected) return false;
    }
    return true;
}