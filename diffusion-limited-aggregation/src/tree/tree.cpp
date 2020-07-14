#include "tree.h"

struct dla_tree* new_dla_tree ()
{
    struct dla_tree *result = (struct dla_tree*)malloc(sizeof(struct dla_tree));

    result->point_list = new_walker_list();
    result->segment_list = new_segment_list();

    return result;    
}

void free_dla_tree (struct dla_tree *the_tree)
{
    if (the_tree->point_list)
        free_walker_list(the_tree->point_list);
    if (the_tree->segment_list)
        free_segment_list(the_tree->segment_list);
    free(the_tree);
}

void grow_tree (struct dla_tree *the_tree, struct user_options *the_options)
{
    struct walker_list *the_point_list = the_tree->point_list;
    struct segment_list *the_segment_list = the_tree->segment_list;

    uint32_t max_number_iterations = the_options->max_num_iter;
    uint32_t seed = the_options->seed;
    double *root_pos = the_options->root_pos;
    bool use_respawn = the_options->use_respawn;
    bool use_initial_network = the_options->use_initial_network;
    char *output_dir = the_options->output_dir;

    create_directory(output_dir);

    // Get the radius of the walkers
    double walker_radius = 6.0;
    get_parameter_value_from_map(the_options->walker_config->params,"walker_radius",&walker_radius);

    srand(seed);

    struct stop_watch solver_time;
    long total_solver_time = 0;
    init_stop_watch(&solver_time);
    start_stop_watch(&solver_time);

//********* MAIN CONFIGURATION BEGIN **********************************************
    // Make the root
    if (use_initial_network)
    {
        make_root_using_initial_network(the_tree,the_options);
    }
    else
    {
        make_root_default(the_tree,the_options);
    }

    // DEBUG
    write_root(the_point_list,the_options->output_dir);

    // Add the Walkers
    struct walker_list *the_others = new_walker_list();
    for (uint32_t i = 0; i < the_options->max_num_walkers; i++)
    {
        struct walker *walker = new_walker(the_options);
        insert_walker_node(the_others,walker);
    }

    // DEBUG
    // Get the reference to the walker draw domain function
    //set_walker_draw_domain_function_fn *draw_domain_function_ptr_2 = the_options->walker_config->draw_domain_function;
    //draw_domain_function_ptr_2(the_options->walker_config,the_options->root_pos);
    //write_list(the_others,0);
    //exit(0);

    // Get the reference to the walker move function
    set_walker_move_function_fn *move_function_ptr = the_options->walker_config->move_function;

    //********* MAIN LOOP STARTS **********************************************
    for (uint32_t iter = 0; iter < max_number_iterations; iter++)
    {
        //printf("Iteration %u of %u\n",iter,max_number_iterations);
        print_progress_bar(iter,max_number_iterations);

        // DEBUG
        //write_list(the_others,iter);

        // Move each Walker using the Random Walk
        struct walker_node *tmp = the_others->list_nodes;
        while (tmp != NULL)
        {

            // Avoid overcrown of segments
            if (the_tree->point_list->num_nodes > MAX_NUMBER_OF_NODES)
            {
                iter = max_number_iterations;
                printf("\n[-] Warning! Maximum number of nodes reached!\n");
                break;
            }
            
            // Move the current walker
            bool has_deleted = false;
            struct walker *cur_walker = tmp->value;
            move_function_ptr(cur_walker,the_options);

            // Check if walker is stuck to a node of the tree
            uint32_t stuck_index = is_stuck(the_tree->point_list,cur_walker);
            if (stuck_index < the_tree->point_list->num_nodes)
            {
                uint32_t new_index = the_tree->point_list->num_nodes;
                struct segment *the_segment = new_segment(stuck_index,new_index);
                insert_segment_node(the_segment_list,the_segment);

                struct walker *the_new_walker = new_walker(cur_walker->pos[0],cur_walker->pos[1],cur_walker->pos[2],cur_walker->radius);
                insert_walker_node(the_point_list,the_new_walker);

                // Delete procedure
                uint32_t delete_index = tmp->id;
                tmp = tmp->next;
                delete_node(the_others,delete_index);
                has_deleted = true;

                order_list(the_others);

                // DEBUG
                //write_current_tree_to_vtk(the_tree,output_dir);
            }

            if (!has_deleted)
            {
                tmp = tmp->next;
            }
                
        }

        // Check if the respawn walker flag is set
        if (use_respawn)
        {
            // Add more walkers as the tree grows ...
            while (the_others->num_nodes < the_options->max_num_walkers)
            {
                struct walker *the_walker = new_walker(the_options,the_others,iter);
                // If it is not possible to add more walker break this loop
                if (the_walker == NULL)
                {
                    //iter = max_number_iterations;
                    //printf("\n[-] Warning! Maximum number of iterations reached for respawn new walker!\n");
                    break;
                }
                insert_walker_node(the_others,the_walker);
            }
        } 
        else
        {
            if (the_others->num_nodes == 0)
            {
                printf("[!] No more walkers left! Leaving main loop!\n");
                break;
            }
        }       
    }
    //********* MAIN LOOP ENDS **********************************************
    total_solver_time = stop_stop_watch(&solver_time);
    printf("\nTotal solver time: %ld μs ==> %ld s\n",total_solver_time,total_solver_time/1000000);

    // Check if the tree is valid
    if (check_dla_tree(the_tree))
    {
        fprintf(stderr,"[-] ERROR! There are non-unique points in the DLA tree!\n");
        exit(EXIT_FAILURE);
    }

    // Print tree information
    print_dla_tree_info(the_tree,output_dir);

    // Write the final network to a file
    write_tree_to_vtk(the_tree,output_dir);

    // Get the reference to the walker draw domain function
    set_walker_draw_domain_function_fn *draw_domain_function_ptr = the_options->walker_config->draw_domain_function;
    draw_domain_function_ptr(the_options->walker_config,the_options->root_pos);

    // Free the walker list
    free_walker_list(the_others);
}

void make_root_default (struct dla_tree *the_tree, struct user_options *the_options)
{
    struct walker_list *the_point_list = the_tree->point_list;
    struct segment_list *the_segment_list = the_tree->segment_list;

    double *root_pos = the_options->root_pos;
    
    // Get the radius of the walkers
    double walker_radius = 6.0;
    get_parameter_value_from_map(the_options->walker_config->params,"walker_radius",&walker_radius);

    // Insert the root point at the position given in the configuration file
    struct walker *root = new_walker(root_pos[0],root_pos[1],root_pos[2],walker_radius);
    insert_walker_node(the_point_list,root);
}

void make_root_using_initial_network (struct dla_tree *the_tree, struct user_options *the_options)
{
    char *initial_network_filename = the_options->initial_network_filename;
    struct walker_list *the_point_list = the_tree->point_list;
    struct segment_list *the_segment_list = the_tree->segment_list;

    double walker_radius = 6.0;
    get_parameter_value_from_map(the_options->walker_config->params,"walker_radius",&walker_radius);

    std::vector<struct network_point> points;
    std::vector<struct network_line> lines;

    // Read the points and lines from the given initial network
    read_initial_network_file(initial_network_filename,points,lines);

    // Insert the points and build a mapping between the initial network and the DLA tree
    uint32_t counter_points = 0;
    std::vector<uint32_t> mapping;
    mapping.assign(points.size(),0);
    for (uint32_t i = 0; i < lines.size(); i++)
    {
        uint32_t src_id = lines[i].src;
        uint32_t dest_id = lines[i].dest;

        struct network_point src = points[src_id];
        struct network_point dest = points[dest_id];

        double u[3], norm;
        calculate_unitary_vector(src.x,src.y,src.z,dest.x,dest.y,dest.z,\
                                u,norm);
        
        uint32_t cur_index = mapping[src_id];
        uint32_t num_walkers_in_line = norm / walker_radius;
        for (uint32_t k = 0; k < num_walkers_in_line; k++)
        {
            double pos[3];
            pos[0] = src.x + k*u[0]*walker_radius;
            pos[1] = src.y + k*u[1]*walker_radius;
            pos[2] = src.z + k*u[2]*walker_radius;

            struct walker *w = new_walker(pos[0],pos[1],pos[2],walker_radius);
            w->stuck_counter = 1;
            insert_walker_node(the_point_list,w);
            counter_points++;

        }
        mapping[dest_id] = counter_points-1;
    }
    // Insert the lines as segments
    counter_points = 0;
    for (uint32_t i = 0; i < lines.size(); i++)
    {
        uint32_t src_id = lines[i].src;
        uint32_t dest_id = lines[i].dest;

        struct network_point src = points[src_id];
        struct network_point dest = points[dest_id];

        double u[3], norm;
        calculate_unitary_vector(src.x,src.y,src.z,dest.x,dest.y,dest.z,\
                                u,norm);
        
        uint32_t cur_index = mapping[src_id];
        uint32_t num_walkers_in_line = norm / walker_radius;
        for (uint32_t k = 0; k < num_walkers_in_line; k++)
        {
            struct segment *the_segment = new_segment(cur_index,counter_points);
            insert_segment_node(the_segment_list,the_segment);
            cur_index = counter_points;
            counter_points++;
        }
    }
}

void build_graph_from_lists (struct walker_list *p_list, struct segment_list *s_list, std::vector< std::vector<uint32_t> > &graph)
{
    // Initialize list of nodes
    uint32_t num_nodes = p_list->num_nodes;
    graph.assign(num_nodes,std::vector<uint32_t>());
    
    // Insert each segment as an edge
    struct segment_node *tmp = s_list->list_nodes;
    while (tmp != NULL)
    {
        uint32_t src_index = tmp->value->src;
        uint32_t dest_index = tmp->value->dest;

        graph[src_index].push_back(dest_index);

        tmp = tmp->next;
    }
}

void print_dla_tree (struct dla_tree *the_tree)
{
    printf("[tree] Point list ...\n");
    print_list(the_tree->point_list);

    printf("[tree] Segment list ...\n");
    print_list(the_tree->segment_list); 
}

void print_dla_tree_info (struct dla_tree *the_tree, const char output_dir[])
{
    struct walker_list *p_list = the_tree->point_list;
    struct segment_list *s_list = the_tree->segment_list;

    // Using the point and segment lists build the corresponding graph
    std::vector< std::vector<uint32_t> > graph;
    build_graph_from_lists(p_list,s_list,graph);

    // Set and open the output information files
    char filename[200];
    sprintf(filename,"%s/segments_length.dat",output_dir);
    FILE *segments_file = fopen(filename,"w+");

    sprintf(filename,"%s/bifurcations_angle.dat",output_dir);
    FILE *bifurcations_file = fopen(filename,"w+");

    uint32_t num_segments = 0;
    uint32_t biff_counter = 0;
    double mean_segment_length = 0.0;
    double mean_bifurcation_angle = 0.0;

    // Use the graph to pass through each segment(edge)
    for (uint32_t i = 0; i < graph.size(); i++)
    {
        uint32_t src_index = i;

        // For each node pass through its edge list
        for (uint32_t j = 0; j < graph[i].size(); j++)
        {
            uint32_t dest_index = graph[i][j];

            // Segment length calculus
            struct walker_node *src_node = search_walker_node(p_list,src_index);
            struct walker_node *dest_node = search_walker_node(p_list,dest_index);

            double segment_length = calculate_euclidean_norm(src_node->value->pos[0],src_node->value->pos[1],src_node->value->pos[2],\
                                                            dest_node->value->pos[0],dest_node->value->pos[1],dest_node->value->pos[2]);
                                                
            if (segment_length > 0.0)
            {
                fprintf(segments_file,"%g\n",segment_length);

                mean_segment_length += segment_length;
                num_segments++;
            }
            
        }

        // Bifurcation angle calculus
        if (graph[i].size() == 2)       // It is a bifurcation
        {
            uint32_t u_index = graph[i][0];
            uint32_t v_index = graph[i][1];

            struct walker_node *i_node = search_walker_node(p_list,i);
            struct walker_node *u_node = search_walker_node(p_list,u_index);
            struct walker_node *v_node = search_walker_node(p_list,v_index);

            double u[3], v[3], norm;
            calculate_unitary_vector(i_node->value->pos[0],i_node->value->pos[1],i_node->value->pos[2],\
                                    u_node->value->pos[0],u_node->value->pos[1],u_node->value->pos[2],u,norm);
            calculate_unitary_vector(i_node->value->pos[0],i_node->value->pos[1],i_node->value->pos[2],\
                                    v_node->value->pos[0],v_node->value->pos[1],v_node->value->pos[2],v,norm);
            
            double angle = calculate_angle_between_vectors(u,v);
            if (!isnan(angle))
            {
                fprintf(bifurcations_file,"%g\n",angle);

                mean_bifurcation_angle += angle;
                biff_counter++;
            }
        }
    }
    mean_segment_length /= (double)num_segments;
    mean_bifurcation_angle /= (double)biff_counter;

    // Convert the length to {mm}
    //mean_segment_length *= 1000;

    fclose(segments_file);
    fclose(bifurcations_file);

    sprintf(filename,"%s/network_info.dat",output_dir);
    FILE *file_info = fopen(filename,"w+");

    fprintf(file_info,"total_num_segments %u\n",num_segments);
    fprintf(file_info,"mean_segment_length %g\n",mean_segment_length);
    fprintf(file_info,"total_num_bifurcations %u\n",biff_counter);
    fprintf(file_info,"mean_bifurcation_angle %g\n",mean_bifurcation_angle);
    
    fclose(file_info);

}

void write_current_tree_to_vtk (struct dla_tree *the_tree, const char output_dir[])
{
    uint32_t num_points = the_tree->point_list->num_nodes;
    uint32_t num_lines = the_tree->segment_list->num_nodes;

    char filename[50];
    sprintf(filename,"%s/tree_%u.vtk",output_dir,num_points);
    FILE *file = fopen(filename,"w+");

    fprintf(file,"# vtk DataFile Version 3.0\n");
    fprintf(file,"Tree\n");
    fprintf(file,"ASCII\n");
    fprintf(file,"DATASET POLYDATA\n");
    fprintf(file,"POINTS %u float\n",num_points);
    
    // Write the points
    struct walker_node *tmp = the_tree->point_list->list_nodes;
    while (tmp != NULL)
    {
        fprintf(file,"%g %g %g\n",tmp->value->pos[0],\
                                tmp->value->pos[1],\
                                tmp->value->pos[2]);
        tmp = tmp->next;
    }

    fprintf(file,"LINES %u %u\n",num_lines,num_lines*3);

    // Write the lines
    struct segment_node *tmp2 = the_tree->segment_list->list_nodes;
    while (tmp2 != NULL)
    {
        fprintf(file,"2 %u %u\n",tmp2->value->src,tmp2->value->dest);
        tmp2 = tmp2->next;
    }
    
    fclose(file);

}

void write_tree_to_vtk (struct dla_tree *the_tree, const char output_dir[])
{
    uint32_t num_points = the_tree->point_list->num_nodes;
    uint32_t num_lines = the_tree->segment_list->num_nodes;

    char filename[50];
    sprintf(filename,"%s/dla_tree.vtk",output_dir,num_points);
    FILE *file = fopen(filename,"w+");

    fprintf(file,"# vtk DataFile Version 3.0\n");
    fprintf(file,"Tree\n");
    fprintf(file,"ASCII\n");
    fprintf(file,"DATASET POLYDATA\n");
    fprintf(file,"POINTS %u float\n",num_points);
    
    // Write the points
    struct walker_node *tmp = the_tree->point_list->list_nodes;
    while (tmp != NULL)
    {
        fprintf(file,"%g %g %g\n",tmp->value->pos[0],\
                                tmp->value->pos[1],\
                                tmp->value->pos[2]);
        tmp = tmp->next;
    }

    fprintf(file,"LINES %u %u\n",num_lines,num_lines*3);

    // Write the lines
    struct segment_node *tmp2 = the_tree->segment_list->list_nodes;
    while (tmp2 != NULL)
    {
        fprintf(file,"2 %u %u\n",tmp2->value->src,tmp2->value->dest);
        tmp2 = tmp2->next;
    }
    
    fclose(file);

}

void write_root (struct walker_list *l, const char output_dir[])
{
    uint32_t num_points = l->num_nodes;

    char filename[50];
    sprintf(filename,"%s/root.vtk",output_dir);
    FILE *file = fopen(filename,"w+");

    fprintf(file,"# vtk DataFile Version 3.0\n");
    fprintf(file,"Walker\n");
    fprintf(file,"ASCII\n");
    fprintf(file,"DATASET POLYDATA\n");
    fprintf(file,"POINTS %u float\n",num_points);
    
    // Write the points
    struct walker_node *tmp = l->list_nodes;
    while (tmp != NULL)
    {
        fprintf(file,"%g %g %g\n",tmp->value->pos[0],\
                                tmp->value->pos[1],\
                                tmp->value->pos[2]);
        tmp = tmp->next;
    }
    fprintf(file,"VERTICES %lu %lu\n",num_points,num_points*2);    
    for (uint32_t i = 0; i < num_points; i++)
        fprintf(file,"1 %u\n",i);
    
    fclose(file);
}

// UNITARY TEST
bool check_dla_tree (struct dla_tree *the_tree)
{
    bool fail = false;

    struct walker_list *p_list = the_tree->point_list;
    struct segment_list *s_list = the_tree->segment_list;

    struct walker_node *tmp = p_list->list_nodes;
    while (tmp != NULL)
    {
        double *p1 = tmp->value->pos;

        struct walker_node *tmp_2 = p_list->list_nodes;
        while (tmp_2 != NULL)
        {
            double *p2 = tmp_2->value->pos;
            if (tmp->id != tmp_2->id)
            {
                double d = calculate_euclidean_norm(p1[0],p1[1],p1[2],p2[0],p2[1],p2[2]);
                if (d == 0.0)
                {
                    printf("[!] WARNING! Nodes %u and %u are the same! Distance = %g\n",tmp->id,tmp_2->id,d);
                    fail = true;
                }
            }
            
            
            tmp_2 = tmp_2->next;
        }
        tmp = tmp->next;
    }
    return fail;
}