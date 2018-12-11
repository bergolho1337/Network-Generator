#include "../include/lsystem.h"

// Build a Lsystem Purkinje network over a endocardium mesh
Lsystem_Generator::Lsystem_Generator (Lsystem_Config *config)
{
	printf("[Lsystem] Building Lsystem network ...\n");

	this->the_purkinje_network = new Graph();
	this->the_miocardium = new Miocardium();

	this->the_miocardium->read_cloud_points(config->miocardium_filename);
	this->the_miocardium->read_terminal_points(config->terminals_filename);

	this->the_miocardium->set_limits();	

	this->initialize_root_point();
	this->initialize_random_array(config->branch_length);
	
	this->make_root(config->branch_length);
	this->grow_network(config);

}

void Lsystem_Generator::grow_network (Lsystem_Config *config)
{
	unsigned int max_iter = config->max_grow_iterations;

	for (unsigned int k = 0; k < max_iter; k++)
	{	
		unsigned int in_the_queue = growing_nodes.size();
		printf("[Lsystem] Iteration %d! Growing branches ...\n",k+1);
		// Take a growing node out of the queue and apply the growing rule on it
		while (in_the_queue > 0)
		{
			Node *ptr = growing_nodes.front();
			growing_nodes.pop();
			
			grow_branch(ptr,config);
			in_the_queue--; 
		}
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

	d_gra[2] = (-sobel[0][2][2] - 2*sobel[0][2][1] - sobel[0][2][0] - 2*sobel[1][2][2] - 4*sobel[1][2][1] - 2*sobel[1][2][0]) -\
		   (sobel[2][2][2] - 2*sobel[2][2][1] - sobel[2][2][0] + sobel[0][0][2] + 2*sobel[0][0][1] + sobel[0][0][0]) +\
		   (2*sobel[1][0][2] + 4*sobel[1][0][1] + 2*sobel[1][0][0] + sobel[2][0][2] + 2*sobel[2][0][1] + sobel[2][0][0]);

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

}

void Lsystem_Generator::grow_branch (Node *ptr, Lsystem_Config *config)
{
	int branch_type = 0;
	double l_bra = config->branch_length;
	double w1 = config->w_1;
	double cube_size = config->sobel_cube_size;
	double tolerance_tree_collision = config->tolerance_tree_collision;
	double tolerance_miocardium_collision = config->tolerance_miocardium_collision;
	double tolerance_terminal_collision = config->tolerance_terminal_collision;
	
	// Calculate the distance gradient by applying the Sobel filter
	double d_gra[3];
	calculate_gradient(ptr,d_gra,cube_size);

	// TODO: Finish this part ...

	/*	
	// Verifica qual o tipo de ramo
	switch (type)
	{
		// Inicio do ramo, gera dois filhos
		case 0: {
					calculateGrowDirection(g,p,d_gra,1);
					calculateGradient(*g,p,d_gra);
					calculateGrowDirection(g,p,d_gra,-1);
					break;
				}
		// Meio do ramo, gera um filho
		case 1: {
					calculateGrowDirection(g,p,d_gra,0);
					break;
				}
		default: break;
	}
	*/
}

void Lsystem_Generator::make_root (const double l_bra)
{
	Graph *pk = this->the_purkinje_network;
	
	double p1[3], p2[3];
	p1[0] = this->root_point[0]; 
	p1[1] = this->root_point[1]; 
	p1[2] = this->root_point[2];

	p2[0] = this->root_point[0]; 
	p2[1] = this->root_point[1]; 
	p2[2] = this->root_point[2] + l_bra;
	
	pk->insert_node_graph(p1);		
	pk->insert_node_graph(p2);
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

	pk->insert_node_graph(pos);
	pk->insert_edge_graph(1,2);

	growing_nodes.push(pk->get_last_node());
}


Lsystem_Config::Lsystem_Config ()
{

}

// TODO: Put this option in the input configuration file
void Lsystem_Generator::initialize_root_point ()
{
	// Benjamin's root position : (158.5,240.28,54.6504)
	// Scaling to the reduced mesh ...
	this->root_point[0] = 39.625;
	this->root_point[1] = 60.07;
	this->root_point[2] = 13.6626;
}

void Lsystem_Generator::initialize_random_array (const double l_bra)
{
	srand(time(NULL));
	default_random_engine generator;
	normal_distribution<double> distribution(0.0,l_bra*SIGMA_RANDOM_DISTRIBUTION);			
	for (int i = 0; i < MAX_SIZE_RAND_ARRAY; i++)
		this->rand_numbers[i] = distribution(generator);
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
	// TODO: Calculate the terminals of the network
	// Imprimir os escalares dos Nodes
	/*
	fprintf(fileVTK,"POINT_DATA %d\n",g->total_nodes);
	fprintf(fileVTK,"SCALARS vm float 1\n");
	fprintf(fileVTK,"LOOKUP_TABLE default\n");
	ptr = g->listNodes;
	while (ptr != NULL)
	{
		if (!ptr->isTerminal)
			fprintf(fileVTK,"0\n");
		else
			fprintf(fileVTK,"1\n");
		ptr = ptr->next;
	}
	*/

	file.close();
}

void Lsystem_Config::print ()
{
	printf("[Lsystem] max_grow_iterations = %u\n",this->max_grow_iterations);
	printf("[Lsystem] branch_length = %.10lf\n",this->branch_length);
	printf("[Lsystem] w_1 = %.10lf\n",this->w_1);
	printf("[Lsystem] tolerance_tree_collision = %.10lf\n",this->tolerance_tree_collision);
	printf("[Lsystem] tolerance_miocardium_collision = %.10lf\n",this->tolerance_miocardium_collision);
	printf("[Lsystem] tolerance_terminal_collision = %.10lf\n",this->tolerance_terminal_collision);
	printf("[Lsystem] miocardium_filename = %s\n",this->miocardium_filename.c_str());
	printf("[Lsystem] terminals_filename = %s\n",this->terminals_filename.c_str());
	
}

