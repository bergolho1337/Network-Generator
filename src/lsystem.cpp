#include "../include/lsystem.h"

// Default constructor 
Lsystem_Config::Lsystem_Config () { }

// Build a Lsystem Purkinje network over a endocardium mesh
Lsystem_Generator::Lsystem_Generator (Lsystem_Config *config)
{
	printf("[Lsystem] Building Lsystem network ...\n");

	this->the_purkinje_network = new Graph();
	this->the_miocardium = new Miocardium();
	this->cont_segments = 1;

	this->the_miocardium->read_cloud_points(config->miocardium_filename);
	this->the_miocardium->read_terminal_points(config->terminals_filename);

	this->the_miocardium->set_limits();	

	this->initialize_root_point(config->root_point);
	this->initialize_random_array(config->branch_length);
	
	this->make_root(config->branch_length);
	this->grow_network(config);

	printf("%s\n",PRINT_LINE);
	printf("[Lsystem] Purkinje structure was build with sucess!\n");
	printf("[Lsystem] There are %d nodes\n",this->the_purkinje_network->get_total_nodes());
	printf("[Lsystem] There are %d edges\n",this->the_purkinje_network->get_total_edges());
	printf("%s\n",PRINT_LINE);

	this->join_terminals();

	// DEBUG
	//this->count_number_terminals();
	//this->the_purkinje_network->print();

}

void Lsystem_Generator::grow_network (Lsystem_Config *config)
{
	unsigned int max_iter = config->max_grow_iterations;

	for (unsigned int k = 0; k < max_iter; k++)
	{	
		unsigned int in_the_queue = growing_nodes.size();
		printf("[Lsystem] Iteration %d! Growing %d branches ...\n",k+1,in_the_queue);
		// Take a growing node out of the queue and apply the growing rule on it
		while (in_the_queue > 0)
		{
			Node *ptr = growing_nodes.front();
			growing_nodes.pop();
			
			grow_branch(ptr,config,0);
			in_the_queue--; 
		}
		//this->the_purkinje_network->print();
	}
}

double Lsystem_Generator::calculate_convolution (double d_gra[], double sobel[3][3][3])
{
	d_gra[0] = (-sobel[0][0][0] + sobel[0][0][2] - 2*sobel[0][1][0] + 2*sobel[0][1][2] - sobel[0][2][0] + sobel[0][2][2]) +\
			(-2*sobel[1][0][0] + 2*sobel[1][0][2] - 4*sobel[1][1][0] + 4*sobel[1][1][2] - 2*sobel[1][2][0] + 2*sobel[1][2][2]) +\
			(-sobel[2][0][0] + sobel[2][0][2] - 2*sobel[2][1][0] + 2*sobel[2][1][2] - sobel[2][2][0] + sobel[2][2][2]);

	d_gra[1] = (sobel[0][0][2] + 2*sobel[1][0][2] + sobel[2][0][2] - sobel[0][2][2] - 2*sobel[1][2][2] - sobel[2][2][2]) +\
			(2*sobel[0][0][1] + 4*sobel[1][0][1] + 2*sobel[2][0][1] - 2*sobel[0][2][1] - 4*sobel[1][2][1] - 2*sobel[2][2][1]) +\
			(sobel[0][0][0] + 2*sobel[1][0][0] + sobel[2][0][0] - sobel[0][2][0] - 2*sobel[1][2][0] - sobel[2][2][0]);

	d_gra[2] = (-sobel[0][2][2] - 2*sobel[0][2][1] - sobel[0][2][0] - 2*sobel[1][2][2] - 4*sobel[1][2][1] - 2*sobel[1][2][0] - sobel[2][2][2] - 2*sobel[2][2][1] - sobel[2][2][0]) + (sobel[0][0][2] + 2*sobel[0][0][1] + sobel[0][0][0] + 2*sobel[1][0][2] + 4*sobel[1][0][1] + 2*sobel[1][0][0] + sobel[2][0][2] + 2*sobel[2][0][1] + sobel[2][0][0]);


	return sqrt(pow(d_gra[0],2)+pow(d_gra[1],2)+pow(d_gra[2],2));
}

bool Lsystem_Generator::is_inside_cube (const Node *ptr,\
			     const double x, const double y, const double z,\
			     const double width, const double lenght, const double height)
{
	// Verify if the point is inside the voxel, which is given by the mini-cube of size (width,lenght,height)
	return (ptr->x <= x + lenght && ptr->x >= x && \
		ptr->y >= y - width && ptr->y <= y && \
	 	ptr->z >= z - height && ptr->z <= z);
}

bool Lsystem_Generator::check_mini_cube (const double x, const double y, const double z,\
			     const double width, const double lenght, const double height)
{
	Graph *pk = this->the_purkinje_network;
	Node *ptr = pk->get_list_nodes();
	
	// For each node of the graph test if at least one is inside the current mini-cube
	while (ptr != NULL)
	{
		// If we found a node that is inside we break the search
		if (is_inside_cube(ptr,x,y,z,width,lenght,height))
			return true;
		ptr = ptr->next;
	}
	
	return false;
}

double Lsystem_Generator::compute_sobel_filter (const Node *p, double sobel[3][3][3], double d_gra[], const double cube_size)
{
	// Dimensions of the mini-cube that surrounds the growing point
	double width = cube_size; 
	double length = cube_size; 
	double height = cube_size;

	double norm = 0.0;
	while (norm == 0.0)
	{
		// Increase the sizes of the mini-cube in the case the gradient is null
		width++; length++; height++;
		// Set the positions of the first mini-cube to be tested
		double x = p->x - ((width*3.0)/2.0); 
		double y = p->y - ((length*3.0)/2.0); 
		double z = p->z - ((height*3.0)/2.0);
	
		for (int i = 0; i < 3; i++)
		{
			z = p->z + ((height*3.0)/2.0);
			for (int j = 0; j < 3; j++)
			{
				y = p->y + ((length*3.0)/2.0);
				for (int k = 0; k < 3; k++)
				{
					// Avoid the central mini-cube
					if (i == 1 && j == 1 && k == 1) 
						sobel[i][j][k] = 0;
					// But, we check the other ones					
					else
					{
						bool tick = check_mini_cube(x,y,z,width,length,height);
						(tick) ? sobel[i][j][k] = 1 : sobel[i][j][k] = 0;
					}
					y -= length;
				}
				z -= height;	
			}
			x += width;
		}
		
		norm = calculate_convolution(d_gra,sobel);	
	}

	// Return the norm of the distance gradient vector
	return norm;
}

void Lsystem_Generator::calculate_gradient (Node *p, double d_gra[], const double cube_size)
{
	double sobel[3][3][3];

	double norm = compute_sobel_filter(p,sobel,d_gra,cube_size);

	// Inverting the value of the distance gradient vector and calculating its unitary direction
	d_gra[0] = -d_gra[0]/norm;
	d_gra[1] = -d_gra[1]/norm;
	d_gra[2] = -d_gra[2]/norm;	

	//printf("d_gra = (%.10lf,%.10lf,%.10lf)\n",d_gra[0],d_gra[1],d_gra[2]);

}

void Lsystem_Generator::rotate_direction (double d[], const double u[], const double teta)
{
	d[0] = u[0]*(cos(teta) + pow(u[0],2)*(1-cos(teta))) +\
		u[1]*(u[0]*u[1]*(1-cos(teta)) - u[2]*sin(teta)) +\
		u[2]*(u[0]*u[2]*(1-cos(teta)) + u[1]*sin(teta));

	d[1] = u[0]*(u[0]*u[1]*(1-cos(teta)) + u[2]*sin(teta)) +\
		u[1]*(cos(teta) + pow(u[1],2)*(1-cos(teta))) +\
		u[2]*(u[1]*u[2]*(1-cos(teta)) - u[0]*sin(teta));

	d[2] = u[0]*(u[0]*u[2]*(1-cos(teta)) - u[1]*sin(teta)) +\
		u[1]*(u[1]*u[2]*(1-cos(teta)) + u[0]*sin(teta)) +\
		u[2]*(cos(teta) + pow(u[2],2)*(1-cos(teta)));
}	

double Lsystem_Generator::calculate_size_branch (const double l_bra)
{
	int i = rand() % MAX_SIZE_RAND_ARRAY;
	double number = rand_numbers[i];
	//return (l_bra + number);
	return l_bra;
}

void Lsystem_Generator::generate_branch (const Node *gnode, const double d[], const Lsystem_Config *config)
{
	Graph *pk = this->the_purkinje_network;

	double l_bra = config->branch_length;
	double tolerance_tree_collision = config->tolerance_tree_collision;
	double tolerance_miocardium_collision = config->tolerance_miocardium_collision;
	double tolerance_terminal_collision = config->tolerance_terminal_collision;


	// Coordinates of the new node
	double d_new[3];
	d_new[0] = gnode->x + calculate_size_branch(l_bra)*d[0];
	d_new[1] = gnode->y + calculate_size_branch(l_bra)*d[1];
	d_new[2] = gnode->z + calculate_size_branch(l_bra)*d[2];

	//printf("Growing node (%.10lf,%.10lf,%.10lf) --- New node (%.10lf,%.10lf,%.10lf)\n",gnode->x,gnode->y,gnode->z,d_new[0],d_new[1],d_new[2]);	

	// !!!! Verify if the point does NOT inflict any restriction  !!!!
	if (!check_terminals(gnode,d_new[0],d_new[1],d_new[2],tolerance_terminal_collision) && \
	    !check_collision_tree(gnode,d_new[0],d_new[1],d_new[2],tolerance_tree_collision) && \
	    !check_collision_miocardium(gnode,d_new[0],d_new[1],d_new[2],tolerance_miocardium_collision) && \
	    !check_limits(gnode,d_new[0],d_new[1],d_new[2]))
	{
		Node *tmp = pk->insert_node_graph(d_new,gnode);
		if (tmp != NULL)
		{
			pk->insert_edge_graph(gnode->id,tmp->id);
			this->cont_segments++;

			// TODO: Try to repeat this 4 times ...
			if (this->cont_segments <= 1) 
				grow_branch(tmp,config,1);
			// Enqueue the last node
			else
				growing_nodes.push(tmp);
		}
	}
}

void Lsystem_Generator::calculate_grow_direction (const Node *gnode, const double d_gra[], const Lsystem_Config *config, const double theta)
{
	double norm;
	double d[3], d_rot[3];

	double w1 = config->w_1;	

	// The vector d_ori it has already been calculated for the Node (insert_node_graph)
	// The vector d_gra it has already been calculated for the Node (calculate_gradient)

	// Calculate the growing direction of the branch by using the rule
	// Rule: d = d_ori + w1*d_gra

	if (theta > 0.0)
	{
		rotate_direction(d_rot,gnode->d_ori,ANGLE);
		for (int i = 0; i < 3; i++)
			d[i] = gnode->d_ori[i] + w1*d_gra[i];	
	}
	else if (theta < 0.0)
	{
		rotate_direction(d_rot,gnode->d_ori,-ANGLE);
		for (int i = 0; i < 3; i++)
			d[i] = gnode->d_ori[i] - w1*d_gra[i];
	}
	else
	{
		for (int i = 0; i < 3; i++)
			d[i] = gnode->d_ori[i] - w1*d_gra[i];
	}
	// Calculate the unitary vector for the direction
	norm = sqrt(pow(d[0],2) + pow(d[1],2) + pow(d[2],2));
	for (int i = 0; i < 3; i++) d[i] /= norm;

	// Create the branch with the calculated direction
	generate_branch(gnode,d,config);

	// Resetar o contador de segmentos
	this->cont_segments = 1;

}

void Lsystem_Generator::grow_branch (Node *gnode, const Lsystem_Config *config, const int branch_type)
{	

	double cube_size = config->sobel_cube_size;

	// Calculate the distance gradient by applying the Sobel filter over the growing node
	double d_gra[3];
	calculate_gradient(gnode,d_gra,cube_size);

	// Verify the type of branch
	switch (branch_type)
	{
		// Generate two offsprings
		case 0: {
				calculate_grow_direction(gnode,d_gra,config,1);
				calculate_gradient(gnode,d_gra,cube_size);
				calculate_grow_direction(gnode,d_gra,config,-1);
				break;
			}
		// Middle of the branch generate only one offspring
		case 1: {
				calculate_grow_direction(gnode,d_gra,config,0);
				break;
			}
		default: break;
	}
}

void Lsystem_Generator::make_root (const double l_bra)
{
	Graph *pk = this->the_purkinje_network;
	
	Node *prev;
	double p1[3], p2[3];
	p1[0] = this->root_point[0]; 
	p1[1] = this->root_point[1]; 
	p1[2] = this->root_point[2];

	p2[0] = this->root_point[0]; 
	p2[1] = this->root_point[1]; 
	p2[2] = this->root_point[2] + l_bra;
	
	prev = pk->insert_node_graph(p1,NULL);		
	prev = pk->insert_node_graph(p2,prev);
	pk->insert_edge_graph(0,1);
		
	// Link the last node to the nearest point in the miocardium
	link_to_miocardium();
	
}

void Lsystem_Generator::link_to_miocardium ()
{
	Graph *pk = this->the_purkinje_network;
	int n = this->the_miocardium->num_cloud_points;
	Point *cloud_points = this->the_miocardium->cloud_points;
	Node *last_node = this->the_purkinje_network->get_last_node();

	double dist = DBL_MAX;
	int id_most_near = -1;
	// Find the nearest point on the miocardium
	for (int i = 0; i < n; i++)
	{
		double norm = calc_norm(last_node->x,last_node->y,last_node->z,cloud_points[i].x,cloud_points[i].y,cloud_points[i].z);
		if (norm < dist && norm > TOLERANCE_NEAREST)
		{
			id_most_near = i;
			dist = norm;
		}
	}
	// Create a new node and make the link
	double pos[3];
	pos[0] = cloud_points[id_most_near].x;
	pos[1] = cloud_points[id_most_near].y;
	pos[2] = cloud_points[id_most_near].z;

	last_node = pk->insert_node_graph(pos,last_node);
	pk->insert_edge_graph(1,2);

	// Repeat this again to build a segment inside the miocardium
	dist = DBL_MAX;
	id_most_near = -1;
	for (int i = 0; i < n; i++)
	{
		double norm = calc_norm(last_node->x,last_node->y,last_node->z,cloud_points[i].x,cloud_points[i].y,cloud_points[i].z);
		if (norm < dist && norm > TOLERANCE_NEAREST)
		{
			id_most_near = i;
			dist = norm;
		}
	}
	pos[0] = cloud_points[id_most_near].x;
	pos[1] = cloud_points[id_most_near].y;
	pos[2] = cloud_points[id_most_near].z;

	last_node = pk->insert_node_graph(pos,last_node);
	pk->insert_edge_graph(2,3);

	growing_nodes.push(pk->get_last_node());
}

void Lsystem_Generator::initialize_root_point (const double root_point[])
{
	// Benjamin's root position : (158.5,240.28,54.6504)
	// Scaling to the reduced mesh ...
	this->root_point[0] = root_point[0];
	this->root_point[1] = root_point[1];
	this->root_point[2] = root_point[2];
}

void Lsystem_Generator::initialize_random_array (const double l_bra)
{
	srand(time(NULL));
	default_random_engine generator;
	normal_distribution<double> distribution(0.0,l_bra*SIGMA_RANDOM_DISTRIBUTION);			
	for (int i = 0; i < MAX_SIZE_RAND_ARRAY; i++)
		this->rand_numbers[i] = distribution(generator);
}

bool Lsystem_Generator::check_collision_tree (const Node *gnode, const double x, const double y, const double z, const double tolerance)
{
	Graph *pk = this->the_purkinje_network;
	
	Node *ptr = pk->get_list_nodes();
	while (ptr->next != NULL)
	{
		if (calc_norm(ptr->x,ptr->y,ptr->z,x,y,z) < tolerance && ptr->is_terminal == false)
		{
			//printf("Collision tree! Between indexes %d -- %d\n",gnode->id,ptr->id);
			return true;
		}
		ptr = ptr->next;
	}
	return false;
}

bool Lsystem_Generator::check_collision_miocardium (const Node *gnode, const double x, const double y, const double z, const double tolerance)
{
	Graph *pk = this->the_purkinje_network;

	int n = this->the_miocardium->num_cloud_points;
	Point *cloud_points = this->the_miocardium->cloud_points;
	
	// Find the nearest point in the cloud of points of the tissue from the growing Node
	int id = search_most_near(cloud_points,n,x,y,z);
	double dist = calc_norm(cloud_points[id].x,cloud_points[id].y,cloud_points[id].z,x,y,z);

	// Check if the distance is in the tolerance
	if (dist < tolerance)
	{
		//printf("Collision miocardium -- Index = %d -- Dist = %.10lf\n",id,dist);
		// Insert the Node from the miocardium
		double pos[3];
		pos[0] = cloud_points[id].x;
		pos[1] = cloud_points[id].y;
		pos[2] = cloud_points[id].z;
		Node *tmp = pk->insert_node_graph(pos,gnode);
		if (tmp != NULL)
		{
			pk->insert_edge_graph(gnode->id,tmp->id);
			// Enqueue ths node ...
			growing_nodes.push(tmp);
			
			return true;
		}
		else
		{
			return false;
		}
	}

	return false;
}

bool Lsystem_Generator::check_limits (const Node *gnode, const double x, const double y, const double z)
{
	double *max_xyz = this->the_miocardium->max_xyz;
	double *min_xyz = this->the_miocardium->min_xyz;
	
	// TODO: Rethink this condition ...
	if (x >= min_xyz[0]-1 && x >= max_xyz[0]+1 && y >= min_xyz[1]-1 && y >= max_xyz[1]+1 \
	    && z >= min_xyz[2]-1 && z >= max_xyz[2]+1)
	{
		return false;
	}
	else
	{
		//printf("Out of limits\n");
		return true;
	}
}

int Lsystem_Generator::search_most_near (const Point *arr, const unsigned int n, const double x, const double y, const double z)
{
	double dist = DBL_MAX;
	int id_nearest = -1;
	for (unsigned int i = 0; i < n; i++)
	{
		if (!arr[i].taken)
		{
			double norm = calc_norm(arr[i].x,arr[i].y,arr[i].z,x,y,z);
			if (norm < dist)
			{
				dist = norm;
				id_nearest = i;
			}
		}
	}
	return id_nearest;
}

// Check if the growing node is not too close of a terminal point
// True = Fail the test
// False = Pass the test
bool Lsystem_Generator::check_terminals(const Node *gnode, const double x, const double y, const double z, const double tolerance)
{
	// Checar se o anterior for terminal retornar erro
	if (gnode->is_terminal) return true;

	Graph *pk = this->the_purkinje_network;

	unsigned int n = this->the_miocardium->num_terminal_points;
	Point *terminals = this->the_miocardium->terminal_points;

	// Find the nearest terminal from the growing node position
	int id = search_most_near(terminals,n,x,y,z);
	double dist = calc_norm(terminals[id].x,terminals[id].y,terminals[id].z,x,y,z);

	// Check if the distance is less than collision tolerance	
	if (dist < tolerance)
	{
		//printf("Terminal colision! Index = %d\n",id);
		// Mark this terminal as taken
		terminals[id].taken = true;
		
		// Copy the coordinate values to an array
		double pos[3];
		pos[0] = terminals[id].x;
		pos[1] = terminals[id].y;
		pos[2] = terminals[id].z;

		// TODO: Find a way to apply the Sobel filter here too ...
		// Insert the terminal node in the graph
		Node *tmp = pk->insert_node_graph(pos,gnode);
		if (tmp != NULL)
		{
			tmp->is_terminal = true;
			pk->insert_edge_graph(gnode->id,tmp->id);
			// !!! We DO NOT enqueue this node for growing
			return true;
		}
		//
		else
		{
			return false;
		}
	}	

	return false;
	
}

void Lsystem_Generator::join_terminals ()
{
	Graph *pk = this->the_purkinje_network;

	int n = this->the_miocardium->num_terminal_points;
	Point *terminals = this->the_miocardium->terminal_points;

	// Pass through all the terminals and verify if there is one that has not been already taken
	for (int i = 0; i < n; i++)
	{
		if (!terminals[i].taken)
		{
			terminals[i].taken = true;
			double dist = DBL_MAX;
			int id_most_near = -1;
			double x, y, z;
			Node *ptr = pk->get_list_nodes();
			while (ptr != NULL)
			{
				double norm = calc_norm(terminals[i].x,terminals[i].y,terminals[i].z,ptr->x,ptr->y,ptr->z);
				if (norm < dist && ptr->num_edges < 2)
				{
					dist = norm;
					id_most_near = ptr->id;
					x = ptr->x;
					y = ptr->y;
					z = ptr->z;	
				}
				ptr = ptr->next;
			}
			// Get the reference to the nearest Node of a terminal
			Node *prev = pk->search_node(id_most_near);
			
			double pos[3];
			pos[0] = terminals[i].x;
			pos[1] = terminals[i].y;
			pos[2] = terminals[i].z;
			Node *tmp = pk->insert_node_graph(pos,prev);

			if (tmp != NULL)
			{
				pk->insert_edge_graph(prev->id,tmp->id);
				tmp->is_terminal = true;
			}
			/*
			printf("=====================================================================\n");
			printf("Most near Node = %d\n",id_most_near);
			printf("Terminal point = (%.10lf,%.10lf,%.10lf) - Most near Node = (%.10lf,%.10lf,%.10lf)\n", \
				terminals[i].x,terminals[i].y,terminals[i].z,x,y,z);
			printf("Distance = %.10lf\n",dist);
			printf("=====================================================================\n");
			*/
		}
	}
}

void Lsystem_Generator::count_number_terminals ()
{
	int n = this->the_miocardium->num_terminal_points;
	Point *terminals = this->the_miocardium->terminal_points;

	int cont = 0;
	for (int i = 0; i < n; i++)
	{
		if (terminals[i].taken) 
			cont++;
	}
	printf("[!] There are %d terminals\n",cont);
}

void Lsystem_Generator::write_network_to_VTK ()
{
	Graph *pk = this->the_purkinje_network;
	std::string filename("networks/tree_nterm.vtk");
	std::ofstream file(filename.c_str());
	Node *ptr;
	Edge *ptrl;

	file << "# vtk DataFile Version 3.0" << std::endl;
	file << "Lsystem" << std::endl;
	file << "ASCII" << std::endl;
	file << "DATASET POLYDATA" << std::endl;
	
	// Write the nodes
	ptr = pk->get_list_nodes();
	file << "POINTS " << pk->get_total_nodes() << " float" << std::endl;
	while (ptr != NULL)
	{
		file << std::setprecision(10) << ptr->x << " " << ptr->y << " " << ptr->z << std::endl;
		ptr = ptr->next;
	}
	
	// Write the edges
	file << "LINES " << pk->get_total_edges() << " " << pk->get_total_edges()*3 << std::endl;
	ptr = pk->get_list_nodes();
	while (ptr != NULL)
	{
		ptrl = ptr->list_edges;
		while (ptrl != NULL)
		{
			file << "2 " << ptr->id << " " << ptrl->id << std::endl;
			ptrl = ptrl->next;
		}
		ptr = ptr->next;
	}

	// Write the scalars (terminals)
	file << "POINT_DATA " << pk->get_total_nodes() << std::endl;
	file << "SCALARS terminal float 1" << std::endl;
	file << "LOOKUP_TABLE default" << std::endl;
	
	ptr = pk->get_list_nodes();
	while (ptr != NULL)
	{
		if (!ptr->is_terminal)
			file << "0" << std::endl;
		else
			file << "1" << std::endl;
		ptr = ptr->next;
	}

	file.close();
}

void Lsystem_Config::print ()
{
	printf("[Lsystem] root_point = (%.10lf,%.10lf,%.10lf)\n",this->root_point[0],this->root_point[1],this->root_point[2]);
	printf("[Lsystem] max_grow_iterations = %u\n",this->max_grow_iterations);
	printf("[Lsystem] branch_length = %.10lf\n",this->branch_length);
	printf("[Lsystem] w_1 = %.10lf\n",this->w_1);
	printf("[Lsystem] tolerance_tree_collision = %.10lf\n",this->tolerance_tree_collision);
	printf("[Lsystem] tolerance_miocardium_collision = %.10lf\n",this->tolerance_miocardium_collision);
	printf("[Lsystem] tolerance_terminal_collision = %.10lf\n",this->tolerance_terminal_collision);
	printf("[Lsystem] miocardium_filename = %s\n",this->miocardium_filename.c_str());
	printf("[Lsystem] terminals_filename = %s\n",this->terminals_filename.c_str());
	
}

