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
    //set_pmj_data();
    
    //test1(this);
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
    for (uint32_t i = 0; i < this->cloud_points.size(); i++)
        delete this->cloud_points[i];
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

    read_cloud_points();
    read_pmj_locations(options->pmj_config);
    //read_obstacles();

    start_stop_watch(&solver_time);

    grow_tree_using_cloud_points(options);

    long res_time = stop_stop_watch(&solver_time);
	double conv_rate = 1000.0*1000.0*60.0;
    printf("[INFO] Resolution time = %ld μs (%g min)\n",res_time,res_time/conv_rate);

    // Unitary test
    bool sucess = check_bifurcation_rule();
    if (!sucess)
    {
        printf("[-] ERROR! Bifurcation rule failure!\n");
        exit(EXIT_FAILURE);
    }

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
    this->cur_rand_index = 0;
    this->max_lat_error = __DBL_MIN__;
    this->min_max_aprox_lat[0] = __DBL_MAX__; this->min_max_aprox_lat[1] = __DBL_MIN__;
    this->min_max_ref_lat[0] = __DBL_MAX__; this->min_max_ref_lat[1] = __DBL_MIN__;
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
    this->cloud_points_filename = options->cloud_points_filename;

    this->Q_TERM = Q_PERF / this->N_term;
    this->R_PERF = powf(((3.0*V_PERF)/(4.0 * M_PI)), (1.0/3.0));
}

void CCO_Network::set_cost_function ()
{
    if (this->cost_function_name == "minimize_custom_function")
    {
        this->cost_fn = new CustomFunction;
    }
    else if (this->cost_function_name == "minimize_activation_time_function")
    {
        this->cost_fn = new ActivationTimeFunction;
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

void CCO_Network::get_segment_length (std::vector<double> &segments)
{
    for (uint32_t i = 0; i < this->segment_list.size(); i++)
    {
        Segment *cur_segment = this->segment_list[i];

        // Segment length (This value is given in {m})
        double length = cur_segment->length;
        if (length == 0.0)
        {
            printf("[!] WARNING! Segment length is equal zero! Indexes: %u %u\n",cur_segment->src->id,cur_segment->dest->id);
        }

        segments.push_back(length * M_TO_MM);
    }
}

void CCO_Network::get_bifurcation_angles(std::vector<double> &angles)
{
    for (uint32_t i = 0; i < this->segment_list.size(); i++)
    {
        Segment *cur_segment = this->segment_list[i];

        if (cur_segment->left != NULL && cur_segment->right != NULL)
        {
            double u[3], v[3], angle;
            cur_segment->left->calc_unitary_vector(u);
            cur_segment->right->calc_unitary_vector(v);
            angle = calc_angle_between_vectors(u,v);

            angles.push_back(angle);
        }
    }
}

void CCO_Network::read_cloud_points ()
{
    if (!this->using_cloud_points)
    {
        fprintf(stderr,"[-] ERROR! You must supply a cloud of points!\n");
        exit(EXIT_FAILURE); 
    }

    bool sucess = read_points_from_vtk(this->cloud_points_filename.c_str(),this->cloud_points);
    this->cloud_points_connected.assign(this->cloud_points.size(),false);
    if (sucess) printf("[cco] Cloud of points was sucessfully loaded from: '%s'!\n",this->cloud_points_filename.c_str());
}

void CCO_Network::read_pmj_locations (PMJConfig *config)
{
    if (this->using_pmj_location)
    {
        if (config->location_filename.size() == 0)
        {
            this->pmj_data = new PMJ();
        }
        else
        {
            this->pmj_data = new PMJ(config);
        }
    }
    else
    {
        this->pmj_data = new PMJ();
    }
}

void CCO_Network::grow_tree_using_cloud_points (User_Options *options)
{
    CostFunctionConfig *cost_function_config = options->cost_function_config;
    LocalOptimizationConfig *local_opt_config = options->local_opt_config;
    //struct pruning_config *pruning_config = options->pruning_config;

    // Root placement
    bool using_initial_network = options->use_initial_network;
    if (using_initial_network)
        make_root_using_initial_network();
    else
        make_root_using_cloud_points();
    
    // MAIN ITERATION LOOP
    while (this->num_terminals < this->N_term)
    {
        printf("%s\n",PRINT_LINE);
        printf("[cco] Working on terminal number %d\n",this->num_terminals+1);
        fprintf(log_file,"%s\n",PRINT_LINE);
        fprintf(log_file,"[cco] Working on terminal number %d\n",this->num_terminals+1);

        if (!this->using_pmj_location)
        {
            bool sucess = generate_terminal_using_cloud_points(cost_function_config,local_opt_config);
            if (sucess)
            {
                write_to_vtk_iteration();
            }
        }
        else
        {
            bool sucess = generate_terminal_using_cloud_points(cost_function_config,local_opt_config);
            if (sucess)
            {
                write_to_vtk_iteration();
            }

            // Active PMJ
            if (this->num_terminals % this->pmj_data->connection_rate == 0)
            {
                bool sucess = attempt_pmj_connection(cost_function_config,local_opt_config);
                if (sucess)
                {
                    write_to_vtk_iteration();
                }
            }
        }

        printf("%s\n",PRINT_LINE);
        fprintf(log_file,"%s\n",PRINT_LINE);
    }

    if (this->using_pmj_location)
    {
        printf("Number of connected PMJ's: %u/%u\n",this->pmj_data->total_num_connected,this->pmj_data->connected.size());

        // Connect the remaining PMJ's
        for (uint32_t i = 0; i < this->pmj_data->points.size(); i++)
        {
            Point *cur_pmj = this->pmj_data->points[i];

            while (!this->pmj_data->connected[i])
            {
                bool sucess = force_pmj_connection(cost_function_config,local_opt_config,cur_pmj);
                if (sucess)
                {
                    write_to_vtk_iteration();
                }   
            }
        }

        // Connect the remaining inactives terminals
        uint32_t num_remanaing_inactives = (this->N_term+this->pmj_data->points.size()) - (this->num_terminals);
        for (uint32_t i = 0; i < num_remanaing_inactives; i++)
        {
            printf("%s\n",PRINT_LINE);
            printf("[cco] Working on inactive terminal number %d\n",this->num_terminals+1);
            fprintf(log_file,"%s\n",PRINT_LINE);
            fprintf(log_file,"[cco] Working on inactive terminal number %d\n",this->num_terminals+1);

            bool sucess = generate_terminal_using_cloud_points(cost_function_config,local_opt_config);
            if (sucess)
            {
                write_to_vtk_iteration();
            }

            printf("%s\n",PRINT_LINE);
            fprintf(log_file,"%s\n",PRINT_LINE);
        }
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
        uint32_t selected_index = sort_point_from_cloud(B);
        double root_length = euclidean_norm(A->x,A->y,A->z,B->x,B->y,B->z);

        //if (root_length >= d_threashold && !has_intersect_obstacle(x_prox,x_inew,obstacle_faces))
        if (root_length >= d_threashold)
        {
            this->cloud_points_connected[selected_index] = true;
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

void CCO_Network::make_root_using_initial_network ()
{
    
    

}

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
    inew->beta = r_left / r_par;
    iconn->beta = r_right / r_par;

    // CONSTANT RADIUS
    //double radius_ratio = r_left / r_right;
    //inew->beta = radius_ratio;
    //iconn->beta = radius_ratio;

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
        ipar_left->beta = r_left / r_par;
        ipar_right->beta = r_right / r_par;
        
        // CONSTANT RADIUS
        //radius_ratio = r_left / r_right;
        //ipar_left->beta = radius_ratio;
        //ipar_right->beta = radius_ratio;

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
    eliminate_point_from_list(M);
    eliminate_point_from_list(T);
    eliminate_segment_from_list(inew);
    eliminate_segment_from_list(ibiff);

    // Update "ndist" from "iconn" until we reach the root
    update_ndist(iconn,true);

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
        uint32_t selected_index = sort_point_from_cloud(new_term);

        // Check the distance criterion for this point
        point_is_ok = connection_search(new_term,d_threash);

        if (point_is_ok) point_is_ok = check_collisions_and_fill_feasible_segments(new_term,feasible_segments);

        if (point_is_ok)
        {
            // COST FUNCTION
            iconn = cost_fn->eval(this,cost_function_config,local_opt_config,\
                                feasible_segments,\
                                new_term);
            
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

            // GIVE UP
            if (d_threash < 1.0E-10)
            {
                delete new_term;

                return sucess;
            }
                
        }
        // The point is a valid one and we can eliminate it from the cloud
        else
        {   
            printf("[!] Selected point %u -- (%g %g %g)\n",selected_index,this->cloud_points[selected_index]->x,this->cloud_points[selected_index]->y,this->cloud_points[selected_index]->z);
            this->cloud_points_connected[selected_index] = true;
        }
    }

    Segment *inew = this->build_segment(local_opt_config,iconn->id,new_term);
    sucess = true;

    if (check_null_segments())
    {
        sucess = false;
        restore_state_tree(iconn);
        exit(1);
    }
    else
    {
        sucess = true;
    }

    delete new_term;

    return sucess;
}

bool CCO_Network::generate_terminal_using_pmj_locations (CostFunctionConfig *cost_function_config, LocalOptimizationConfig *local_opt_config, Point *pmj_point, const bool evaluate)
{
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

    while (!point_is_ok)
    {
        // Reset the feasible segments list
        feasible_segments.clear();

        // Check the distance criterion for this point
        point_is_ok = connection_search(pmj_point,d_threash);

        if (point_is_ok) point_is_ok = check_collisions_and_fill_feasible_segments(pmj_point,feasible_segments);

        if (point_is_ok)
        {

            // [COST FUNCTION]
            iconn = cost_fn->eval(this,cost_function_config,local_opt_config,\
                                feasible_segments,\
                                pmj_point);
            
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

                counter_tries++;
            }

            // GIVE UP
            if (counter_tries > this->pmj_data->max_connection_tries)
                return sucess;
                
        }
        // The point is a valid one and we can eliminate it from the cloud
        else
        {

        }
    }

    // Build the new segment of the tree
    Segment *inew = this->build_segment(local_opt_config,iconn->id,pmj_point);

    // [EVALUATE PMJ ACTIVATION]
    if (evaluate)
    {
        // Step 1: Try to connect using normal CCO
        sucess = evaluate_pmj_local_activation_time(inew,pmj_point,cost_function_config);
        //if (sucess) printf("[+] PMJ point %u was connected in step 1!\n",pmj_point->id);

        // Step 2: Try to connect using a region radius
        if (!sucess)
        {   
            sucess = attempt_connect_using_region_radius(pmj_point,cost_function_config);
            //if (sucess) printf("[+] PMJ point %u was connected in step 2!\n",pmj_point->id);
        }
    }
    else
    {
        double lat = calc_terminal_local_activation_time(inew) + this->lat_offset;

        // Compute the LAT error for the current PMJ
        double lat_error = (lat - pmj_point->lat);
        this->pmj_data->error[pmj_point->id] = lat_error;
        this->pmj_data->aprox[pmj_point->id] = lat;
        
        // Mark this PMJ as untouchable
        inew->can_touch = false;
        
        sucess = true;
    }
    
    
    return sucess;
}

bool CCO_Network::connection_search (Point *p, const double d_threash)
{
    // Copy the coordinates of the new terminal to an array
    double pos[3];
    pos[0] = p->x;
    pos[1] = p->y;
    pos[2] = p->z;

    for (uint32_t i = 0; i < this->segment_list.size(); i++)
    {
        Segment *s = this->segment_list[i];

        if (!distance_criterion(s,pos,d_threash))
            return false;
    }
    return true;
}

bool CCO_Network::distance_criterion (Segment *s, const double pos[], const double d_threash)
{
    double d_proj = calc_dproj(s,pos);

    double d_crit;
    if (d_proj >= 0 && d_proj <= 1)
        d_crit = calc_dortho(s,pos);
    else
        d_crit = calc_dend(s,pos);

    return (d_crit < d_threash) ? false : true;
}

bool CCO_Network::check_collisions_and_fill_feasible_segments (Point *p, std::vector<Segment*> &feasible_segments)
{
    fprintf(log_file,"[!] Checking collisions!\n");

    std::vector< std::pair<double,uint32_t> > arr;
    for (uint32_t i = 0; i < this->segment_list.size(); i++)
    {
        Segment *s = this->segment_list[i];

        if (!has_collision(s,p))
        {
            double middle_pos[3];
            s->calc_middle_point(middle_pos);

            double dist = euclidean_norm(middle_pos[0],middle_pos[1],middle_pos[2],\
                                    p->x,p->y,p->z);

            arr.push_back(std::make_pair(dist,s->id));

            feasible_segments.push_back(s);
        }      
    }

    // Sort teh segments by the distance to the middle point
    std::sort(arr.begin(),arr.end());

    // Get only the closest NCONN segments
    for (uint32_t i = 0; i < arr.size() && feasible_segments.size() < NCONN; i++)
    {
        uint32_t id = arr[i].second;
        
        Segment *s = this->segment_list[id];

        feasible_segments.push_back(s);
    }

    return (feasible_segments.size() == 0) ? false : true;
}

bool CCO_Network::check_bifurcation_rule ()
{
    double gamma = this->gamma;
    for (uint32_t i = 0; i < this->segment_list.size(); i++)
    {
        Segment *cur_segment = this->segment_list[i];
        Segment *left = cur_segment->left;
        Segment *right = cur_segment->right;

        if (left != NULL && right != NULL)
        {
            double r = pow(cur_segment->radius,gamma);
            double r_left = pow(left->radius,gamma);
            double r_right = pow(right->radius,gamma);

            double diff = r_left + r_right - r;
            fprintf(log_file,"%g = %g + %g --> %g = %g\n",r,r_left,r_right,r,r_left + r_right);

            if (diff > 1.0E-8)
            {
                printf("[-] ERROR! Bifurcation rule was not valid !\n");
                return false;
            }
        }
    }
    return true;
}

bool CCO_Network::check_null_segments ()
{
    bool result = false;
    for (uint32_t i = 0; i < this->segment_list.size(); i++)
    {
        if (this->segment_list[i]->length == 0.0)
        {
            printf("[!] WARNING! Segment '%u' has a null length!\n",this->segment_list[i]->id);
            result = true;
        } 
    }
    return result;
}

// p1 -> proximal || p2 -> distal || p3 -> middle || p4 -> new 
bool CCO_Network::has_collision (Segment *s, Point *p)
{
    double middle_pos[3];
    s->calc_middle_point(middle_pos);

    //printf("[+] Trying connection with segment %u\n",s->id);

    for (uint32_t i = 0; i < this->segment_list.size(); i++)
    {
        Segment *cur_segment = this->segment_list[i];
        uint32_t cur_id = cur_segment->id;

        if (s->id != cur_id)
        {
            //printf("\t[!] Checking collison between segment %d\n",cur_id);

            Point *src = cur_segment->src;
            Point *dest = cur_segment->dest;

            bool intersect = collision_detection(src->x,src->y,src->z,\
                                            dest->x,dest->y,dest->z,\
                                            middle_pos[0],middle_pos[1],middle_pos[2],\
                                            p->x,p->y,p->z);
            if (intersect)
            {
                printf("\t[-] ERROR! Intersection with segment %d !\n",cur_id);
                return true;
            }
        }
    }

    return false;
}

bool CCO_Network::evaluate_pmj_local_activation_time (Segment *inew, Point *pmj_point, CostFunctionConfig *cost_function_config)
{
    bool sucess = false;

    // Calculate the LAT for the new PMJ point included in the network
    double lat = calc_terminal_local_activation_time(inew) + this->lat_offset;

    // Compute the LAT error for the current PMJ
    double lat_error = (lat - pmj_point->lat);

    // Rodrigo's suggestion
    //if (error > LAT_ERROR_LIMIT)

    if (fabs(lat_error) < this->pmj_data->lat_error_tolerance)
    {
        sucess = true;
        inew->can_touch = false;
        this->pmj_data->error[pmj_point->id] = lat_error;
        this->pmj_data->aprox[pmj_point->id] = lat;
        //write_pathway(pmj_point);
    }
    else
    {
        prune_segment(inew);
    }

    return sucess;
}

Point* CCO_Network::generate_bifurcation_node (Segment *iconn, LocalOptimizationConfig *local_opt_config)
{
    uint32_t cur_num_nodes = this->point_list.size(); 
    Point *result = new Point(cur_num_nodes);

    if (this->using_local_optimization)
    {
        // The local optimization was not executed yet
        if (local_opt_config->first_call)
        {
            double middle_pos[3];
            iconn->calc_middle_point(middle_pos);
            result->setCoordinate(middle_pos);
        }
        // The best bifurcation position is already stored in the 'best_pos' array
        else
        {
            double *best_pos = local_opt_config->best_pos;
            result->setCoordinate(best_pos);

            // Reset the 'first_call' flag for the next point
            local_opt_config->first_call = true;
        }
    }
    else
    {
        double middle_pos[3];
        iconn->calc_middle_point(middle_pos);
        result->setCoordinate(middle_pos);
    }
    result->setLAT(iconn->calc_middle_point_lat());
    result->setActive(false);

    return result;
}

Point* CCO_Network::generate_terminal_node (Point *p)
{
    uint32_t cur_num_nodes = this->point_list.size();
    Point *result = new Point(cur_num_nodes);

    double pos[3];
    pos[0] = p->x;
    pos[1] = p->y;
    pos[2] = p->z;

    result->setCoordinate(pos);
    result->setActive(p->is_active);
    result->setLAT(p->lat);

    return result;
}

void CCO_Network::update_segment_pointers (Segment *iconn, Segment *ibiff, Segment *inew, const bool is_restore)
{
    // Update pointers when building a new segment
    if (!is_restore)
    {
        // iconn:
        iconn->parent = ibiff;
        iconn->src = inew->src;

        // ibiff:
        ibiff->left = inew;     // CONVENTION: Left will always point to terminal
        ibiff->right = iconn;   // CONVENTION: Right will always point to subtree
        if (ibiff->parent != NULL)
        {
            Segment *ibiff_parent = ibiff->parent;
            if (ibiff_parent->left->id == iconn->id)
                ibiff_parent->left = ibiff;
            if (ibiff_parent->right->id == iconn->id)
                ibiff_parent->right = ibiff;
        }
    }
    // Update pointers when restoring the tree state
    else
    {
        Segment *ibiff_par = ibiff->parent;

        // iconn:
        iconn->parent = ibiff->parent;
        iconn->src = ibiff->src;
        iconn->beta = ibiff->beta;

        // ibiff_par:
        if (ibiff_par != NULL)
        {
            if (ibiff_par->right == ibiff)
                ibiff_par->right = iconn;
            if (ibiff_par->left == ibiff)
                ibiff_par->left = iconn;
        }
    }
}

void CCO_Network::update_ndist (Segment *s, const bool is_restore)
{
    // Update pointers when building a new segment
    if (!is_restore)
    {
        Segment *tmp = s;
        while (tmp->parent != NULL)
        {
            tmp->ndist++;
            tmp = tmp->parent;
        }
        tmp->ndist++;

        // Update the current number of terminals in the network
        this->num_terminals = tmp->ndist;
    }
    // Update pointers when restoring the tree state
    else
    {
        Segment *tmp = s->parent;
        if (tmp != NULL)
        {
            while (tmp->parent != NULL)
            {
                tmp->ndist--;
                tmp = tmp->parent;
            }
            tmp->ndist--;
            this->num_terminals = tmp->ndist;
        }
        else
        {
            this->num_terminals = s->ndist;
        }
    }
}

Segment* CCO_Network::build_segment (LocalOptimizationConfig *local_opt_config, const uint32_t index, Point *new_term)
{
    uint32_t cur_num_segments = this->segment_list.size();
    Segment *iconn = this->segment_list[index];
    
    // Create the bifurcation point
    Point *M = generate_bifurcation_node(iconn,local_opt_config);
    this->point_list.push_back(M);

    // Create the terminal point
    Point *T = generate_terminal_node(new_term);
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
    update_ndist(ibiff,false);

    // Rescale the network from the segments of the new bifurcation
    rescale_tree(ibiff,iconn,inew);

    // Recalculate the segment radius using the 'beta' values
    recalculate_radius();
    recalculate_length();

    // Return a reference to the new node
    return inew;
}

double CCO_Network::calc_terminal_local_activation_time (Segment *term)
{
    // TODO: Consider the case where we don't have a constant valocity
    const double cv = 1900.0; // Conduction velocity across the Purkinje fiber {um/ms}

    //double delta_s = term->length;   // During runtime the network length is given in {m}
    double delta_s = 0.0;   // During runtime the network length is given in {m}
    Segment *tmp = term;
    while (tmp != NULL)
    {
        delta_s += tmp->length;

        tmp = tmp->parent;
    }

    // Return the LAT in {ms}
    return delta_s*M_TO_UM/cv; 
}

uint32_t CCO_Network::sort_point_from_cloud (Point *p)
{
    uint32_t offset = rand() % this->max_rand_offset + 1;

    // Reset the counter
    if (this->cur_rand_index > this->cloud_points.size()-1) 
        this->cur_rand_index = this->cur_rand_index % (this->cloud_points.size()-1);

    uint32_t selected_index = this->cur_rand_index;
    
    // Get the current cloud point data
    bool is_active;
    double pos[3], ref_lat;
    pos[0] = this->cloud_points[selected_index]->x;
    pos[1] = this->cloud_points[selected_index]->y;
    pos[2] = this->cloud_points[selected_index]->z;
    ref_lat = this->cloud_points[selected_index]->lat;
    is_active = this->cloud_points[selected_index]->is_active;

    // Set the point data
    p->setCoordinate(pos);
    p->setLAT(ref_lat);
    p->setActive(is_active);
    
    // Increase the counter for the next cloud point
    this->cur_rand_index += offset;

    return selected_index;
}

void CCO_Network::eliminate_point_from_list (Point *p)
{
    this->point_list.erase(this->point_list.begin() + p->id);
    order_point_list();
    delete p;
    p = NULL;
}

void CCO_Network::eliminate_segment_from_list (Segment *s)
{
    this->segment_list.erase(this->segment_list.begin() + s->id);
    order_segment_list();
    delete s;
    s = NULL;
}

void CCO_Network::order_point_list ()
{
    uint32_t counter = 0;
    for (uint32_t i = 0; i < this->point_list.size(); i++)
    {
        this->point_list[i]->id = counter;
        counter++;
    }
}

void CCO_Network::order_segment_list ()
{
    uint32_t counter = 0;
    for (uint32_t i = 0; i < this->segment_list.size(); i++)
    {
        this->segment_list[i]->id = counter;
        counter++;
    }
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
    eliminate_point_from_list(M);
    eliminate_point_from_list(T);
    eliminate_segment_from_list(inew);
    eliminate_segment_from_list(ibiff);

    // Update "ndist" from "iconn" until we reach the root
    update_ndist(iconn,true);

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
    get_segment_length(segments);

    std::vector<double> angles;
    get_bifurcation_angles(angles);

    double mean_segment_length, std_segment_length;
    calc_mean_std(segments,mean_segment_length,std_segment_length);
    write_vector_to_file(segments,this->output_dir + "/segments_length.dat");

    double mean_biff_angle, std_biff_angle;
    calc_mean_std(angles,mean_biff_angle,std_biff_angle);
    write_vector_to_file(angles,this->output_dir + "/bifurcation_angle.dat");

    write_info_to_file(this->output_dir + "/network_info.txt",\
                    segments,mean_segment_length,std_segment_length,\
                    angles,mean_biff_angle,std_biff_angle);
    
    printf("[INFO] Total number of segment = %u\n",segments.size());
    printf("[INFO] Segment length = %g +/- %g mm\n",mean_segment_length,std_segment_length);
    printf("[INFO] Total number of bifurcations = %u\n",angles.size());
    printf("[INFO] Bifurcation angle = %g +/- %g degrees\n",mean_biff_angle,std_biff_angle);
    if (this->using_pmj_location) 
    {
        printf("[INFO] Number of PMJ's connected = %u/%u\n",this->pmj_data->total_num_connected,this->pmj_data->points.size());
        write_vector_to_file(this->pmj_data->error,this->output_dir + "/pmj_error.dat");

        get_electric_error();
        printf("[INFO] Max Error = %g || Ref Min LAT = %g || Ref Max LAT = %g || Aprox Min LAT = %g || Aprox Max LAT = %g\n",this->max_lat_error,this->min_max_ref_lat[0],this->min_max_ref_lat[1],\
                                                                                                                    this->min_max_aprox_lat[0],this->min_max_aprox_lat[1]);
    }
}

bool CCO_Network::attempt_pmj_connection (CostFunctionConfig *cost_function_config, LocalOptimizationConfig *local_opt_config)
{
    uint32_t cur_num_connected_pmjs = 0;
    uint32_t prev_num_connected_pmjs = 0;
    uint32_t total_num_pmjs_connected_before = this->pmj_data->total_num_connected;
    do
    {
        prev_num_connected_pmjs = cur_num_connected_pmjs;
        for (uint32_t i = 0; i < this->pmj_data->points.size(); i++)
        {
            printf("[cco] Trying to insert PMJ point %u ...\n",i+1);
            fprintf(log_file,"[cco] Trying to insert PMJ point %u ...\n",i+1);

            Point *cur_pmj_point = this->pmj_data->points[i];

            if (!this->pmj_data->connected[i])
            {
                bool sucess = generate_terminal_using_pmj_locations(cost_function_config,local_opt_config,cur_pmj_point,true);
                this->pmj_data->connected[i] = sucess;
                
                if (sucess)
                {
                    printf("\t[!] SUCESS!\n");
                    cur_num_connected_pmjs++;
                    this->pmj_data->total_num_connected++;
                    write_to_vtk_iteration();
                }
            }
            //printf("\t[cco] Previous connected PMJ = %u || Current connected PMJ = %u\n",prev_num_connected_pmjs,cur_num_connected_pmjs);
            //fprintf(log_file,"\t[cco] Previous connected PMJ = %u || Current connected PMJ = %u\n",prev_num_connected_pmjs,cur_num_connected_pmjs);
        }
    } while (prev_num_connected_pmjs != cur_num_connected_pmjs);

    //printf("[cco] Number of connected PMJ in this pass: %u\n",this->total_num_pmjs_connected-total_num_pmjs_connected_before);
    //fprintf(log_file,"[cco] Number of connected PMJ in this pass: %u\n",this->total_num_pmjs_connected-total_num_pmjs_connected_before);
    
    return (this->pmj_data->total_num_connected > total_num_pmjs_connected_before) ? true : false;
}

bool CCO_Network::attempt_connect_using_region_radius (Point *pmj_point, CostFunctionConfig *cost_function_config)
{
    bool sucess = false;
    double center[3], ori_pos[3];
    double radius = this->pmj_data->region_radius;

    center[0] = pmj_point->x;
    center[1] = pmj_point->y;
    center[2] = pmj_point->z;

    for (uint32_t i = 0; i < this->segment_list.size(); i++)
    {
        Segment *cur_segment = this->segment_list[i];
        if (cur_segment->is_terminal() && cur_segment->is_inside_region(center,radius) && cur_segment->is_touchable() )
        {
            // Try to connect the PMJ location directly to a segment by changing the destination node coordinates
            //printf("Segment %u is inside region for PMJ %u\n",cur_segment->id,pmj_point->id);

            // Save the original position
            ori_pos[0] = cur_segment->dest->x;
            ori_pos[1] = cur_segment->dest->y;
            ori_pos[2] = cur_segment->dest->z;

            // Change coordinate position
            cur_segment->dest->x = pmj_point->x;
            cur_segment->dest->y = pmj_point->y;
            cur_segment->dest->z = pmj_point->z;
            cur_segment->length = cur_segment->calc_length();
            cur_segment->can_touch = false;

            double lat = calc_terminal_local_activation_time(cur_segment) + this->lat_offset;
            double lat_error = (lat - pmj_point->lat);

            if (fabs(lat_error) < this->pmj_data->lat_error_tolerance)
            {
                sucess = true;
                this->pmj_data->error[pmj_point->id] = lat_error;
                this->pmj_data->aprox[pmj_point->id] = lat;
                //print_pathway(tmp,pmj_id);

                return sucess;
            }
            else
            {
                // Restore the state
                cur_segment->dest->x = ori_pos[0];
                cur_segment->dest->y = ori_pos[1];
                cur_segment->dest->z = ori_pos[2];
                cur_segment->length = cur_segment->calc_length();
                cur_segment->can_touch = true;
            }
        }
    }
    return sucess;
}

bool CCO_Network::force_pmj_connection (CostFunctionConfig *cost_function_config, LocalOptimizationConfig *local_opt_config, Point *pmj_point)
{
    bool sucess = false;
    double original_value = this->pmj_data->lat_error_tolerance;

    printf("[cco] Forcing connection for PMJ point %u ...\n",pmj_point->id);
    fprintf(log_file,"[cco] Forcing connection for PMJ point %u ... ...\n",pmj_point->id);

    while (!sucess)
    {
        sucess = generate_terminal_using_pmj_locations(cost_function_config,local_opt_config,pmj_point,true);
        if (sucess)
        {
            this->pmj_data->connected[pmj_point->id] = sucess;
            this->pmj_data->total_num_connected++;
        }
        else
        {
            printf("\t[cco] Looseing LAT tolerance error from %g to %g\n",this->pmj_data->lat_error_tolerance,this->pmj_data->lat_error_tolerance*1.2);
            this->pmj_data->lat_error_tolerance *= 1.2;
        }

        // Brute force connection
        if (this->pmj_data->lat_error_tolerance > PMJ_LOOSE_THREASHOLD)
        {
            printf("\t[cco] Forcing connection!\n");
            sucess = generate_terminal_using_pmj_locations(cost_function_config,local_opt_config,pmj_point,false);
            if (!sucess)
            {
                sucess = force_connection_to_closest_segment(cost_function_config,local_opt_config,pmj_point);
            }
            this->pmj_data->connected[pmj_point->id] = sucess;
            if (sucess)
            {
                this->pmj_data->total_num_connected++;
            }
        }
    }

    // Reset to the original value
    this->pmj_data->lat_error_tolerance = original_value;

    return sucess;
}

bool CCO_Network::force_connection_to_closest_segment (CostFunctionConfig *cost_function_config, LocalOptimizationConfig *local_opt_config, Point *pmj_point)
{
    Segment *closest = NULL;
    double min_dist = __DBL_MAX__;

    for (uint32_t i = 0; i < this->segment_list.size(); i++)
    {
        Segment *cur_segment = this->segment_list[i];

        if (cur_segment->is_terminal() && cur_segment->is_touchable())
        {
            double dist = euclidean_norm(cur_segment->dest->x,cur_segment->dest->y,cur_segment->dest->z,\
                                        pmj_point->x,pmj_point->y,pmj_point->z);
            if (dist < min_dist)
            {
                min_dist = dist;
                closest = cur_segment;
            }
        }
    }

    // Change the coordinates
    closest->dest->x = pmj_point->x;
    closest->dest->y = pmj_point->y;
    closest->dest->z = pmj_point->z;
    closest->length = closest->calc_length();
    closest->can_touch = false;
    
    // Calculate the LAT error
    double lat = calc_terminal_local_activation_time(closest) + this->lat_offset;
    double lat_error = fabs(lat-pmj_point->lat);
    this->pmj_data->error[pmj_point->id] = lat_error;
    this->pmj_data->aprox[pmj_point->id] = lat;
    
    return true;
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

Segment* CCO_Network::get_terminal (const double pos[])
{
    uint32_t min_id = 0;
    double min_dist = __DBL_MAX__;

    // Find the index of the closest point to the one given
    for (uint32_t i = 0; i < this->point_list.size(); i++)
    {
        double middle_pos[3];
        middle_pos[0] = this->point_list[i]->x*M_TO_UM;
        middle_pos[1] = this->point_list[i]->y*M_TO_UM;
        middle_pos[2] = this->point_list[i]->z*M_TO_UM;

        double dist = euclidean_norm(middle_pos[0],middle_pos[1],middle_pos[2],\
                                    pos[0],pos[1],pos[2]);
        
        if (dist < min_dist)
        {
            min_dist = dist;
            min_id = i;
        }
    }

    for (uint32_t i = 0; i < this->segment_list.size(); i++)
    {
        if (this->segment_list[i]->dest->id == min_id)
            return this->segment_list[i];
    }
    for (uint32_t i = 0; i < this->segment_list.size(); i++)
    {
        if (this->segment_list[i]->src->id == min_id)
            return this->segment_list[i];
    }
    return NULL;
}

Point* CCO_Network::search_point (const uint32_t index)
{
    return this->point_list[index];
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
    result->max_lat_error = __DBL_MIN__;
    result->min_max_aprox_lat[0] = __DBL_MAX__; result->min_max_aprox_lat[1] = __DBL_MIN__;
    result->min_max_ref_lat[0] = __DBL_MAX__; result->min_max_ref_lat[1] = __DBL_MIN__;
    srand(result->seed);

    memcpy(result->root_pos,this->root_pos,sizeof(double)*3);

    result->using_only_murray_law = this->using_only_murray_law;
    result->start_radius = this->start_radius;
    
    result->using_cloud_points = this->using_cloud_points;
    result->using_local_optimization = this->using_local_optimization;
    result->using_pmj_location = this->using_pmj_location;

    if (this->pmj_data)
    {
        result->pmj_data = this->pmj_data->copy();
    }
        

    memcpy(result->min_max_aprox_lat,this->min_max_aprox_lat,sizeof(double)*2);
    memcpy(result->min_max_ref_lat,this->min_max_ref_lat,sizeof(double)*2);
    
    result->output_dir = "outputs/copy";
    result->cost_function_name = this->cost_function_name;
    result->local_optimization_function_name = this->local_optimization_function_name;
    result->cloud_points_filename = this->cloud_points_filename;

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
    for (uint32_t i = 0; i < this->segment_list.size(); i++)
    {
        Segment *left = this->segment_list[i]->left;
        Segment *right = this->segment_list[i]->right;
        Segment *parent = this->segment_list[i]->parent;

        result->segment_list[i]->left = (left) ? result->segment_list[left->id] : NULL;
        result->segment_list[i]->right = (right) ? result->segment_list[right->id] : NULL;
        result->segment_list[i]->parent = (parent) ? result->segment_list[parent->id] : NULL;
    }

    // Copy cloud of points
    for (uint32_t i = 0; i < this->cloud_points.size(); i++)
    {
        Point *p = new Point(this->cloud_points[i]);
        result->cloud_points.push_back(p);
    }
    
    return result;
}

CCO_Network* CCO_Network::concatenate (CCO_Network *input)
{
    uint32_t offset_points, offset_segments;

    CCO_Network *result = this->copy();

    // Sum the number of terminals from both networks
    result->num_terminals += input->num_terminals;
    result->using_pmj_location |= input->using_pmj_location;

    // Adjust the Min/Max PMJ errors
    result->min_max_aprox_lat[0] = (input->min_max_aprox_lat[0] < this->min_max_aprox_lat[0]) ? input->min_max_aprox_lat[0] : this->min_max_aprox_lat[0];
    result->min_max_aprox_lat[1] = (input->min_max_aprox_lat[1] > this->min_max_aprox_lat[1]) ? input->min_max_aprox_lat[1] : this->min_max_aprox_lat[1];
    result->min_max_ref_lat[0] = (input->min_max_ref_lat[0] < this->min_max_ref_lat[0]) ? input->min_max_ref_lat[0] : this->min_max_ref_lat[0];
    result->min_max_ref_lat[1] = (input->min_max_ref_lat[1] > this->min_max_ref_lat[1]) ? input->min_max_ref_lat[1] : this->min_max_ref_lat[1];
    result->max_lat_error = (input->max_lat_error > this->max_lat_error) ? input->max_lat_error : this->max_lat_error;

    // Concatenate points
    offset_points = result->point_list.size();
    for (uint32_t i = 0; i < input->point_list.size(); i++)
    {
        Point *p = new Point(input->point_list[i]);
        p->id += offset_points;
        result->point_list.push_back(p);
    }

    // Concatenate segments
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
    offset_points = result->cloud_points.size();
    for (uint32_t i = 0; i < input->cloud_points.size(); i++)
    {
        Point *p = new Point(input->cloud_points[i]);
        p->id += offset_points;
        result->cloud_points.push_back(p);
    }

    // Concatenate PMJ's
    if (input->pmj_data)
    {
        result->pmj_data->total_num_connected += input->pmj_data->total_num_connected;
        
        offset_points = result->pmj_data->points.size();
        for (uint32_t i = 0; i < input->pmj_data->points.size(); i++)
        {
            Point *p = new Point(input->pmj_data->points[i]);
            bool connected = input->pmj_data->connected[i];
            double aprox = input->pmj_data->aprox[i];
            double error = input->pmj_data->error[i];
            p->id += offset_points;

            result->pmj_data->points.push_back(p);
            result->pmj_data->connected.push_back(connected);
            result->pmj_data->aprox.push_back(aprox);
            result->pmj_data->error.push_back(error);
        }
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

void CCO_Network::get_electric_error ()
{
    for (uint32_t i = 0; i < this->pmj_data->aprox.size(); i++)
    {
        double lat = this->pmj_data->aprox[i];
        if (lat < this->min_max_aprox_lat[0]) this->min_max_aprox_lat[0] = lat;
        if (lat > this->min_max_aprox_lat[1]) this->min_max_aprox_lat[1] = lat;
    }

    for (uint32_t i = 0; i < this->pmj_data->points.size(); i++)
    {
        double lat = this->pmj_data->points[i]->lat;
        if (lat < this->min_max_ref_lat[0]) this->min_max_ref_lat[0] = lat;
        if (lat > this->min_max_ref_lat[1]) this->min_max_ref_lat[1] = lat;
    }

    for (uint32_t i = 0; i < this->pmj_data->error.size(); i++)
    {
        double error = fabs(this->pmj_data->error[i]);

        if (error > this->max_lat_error) this->max_lat_error = error;
    }
}