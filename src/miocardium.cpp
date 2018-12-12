#include "../include/miocardium.h"

// TODO: Pur this in a separate file 'miocardium.cpp'
Miocardium::Miocardium ()
{
	this->num_cloud_points = 0;
	this->num_terminal_points = 0;
	this->cloud_points = NULL;
	this->terminal_points = NULL;
}

Point::Point () 
{ 
	this->taken = false;
}

void Miocardium::read_cloud_points (const std::string miocardium_filename)
{
	printf("[Miocardium] Reading cloud of points ...\n");

	std::ifstream file;
	std::string str;
	int n;
	double x, y, z, scalar;
	double center_x, center_y, center_z;
	double rootPoint[3];

	file.open(miocardium_filename.c_str());
	if (!file) 
	{
		printf("[-] ERROR! Cannot load miocardium file '%s'!\n",miocardium_filename.c_str());
		exit(EXIT_FAILURE);
	}

	center_x = 0.0;
	center_y = 0.0;
	center_z = 0.0;

	// Read points
	while (file >> str && str != "POINTS");

	file >> n >> str;
	this->num_cloud_points = n;
	this->cloud_points = new Point[n];

	Point *cloud_points = this->cloud_points;
	
	for (int i = 0; i < n; i++)
	{
		file >> x >> y >> z;

		cloud_points[i].x = x;		
		cloud_points[i].y = y;
		cloud_points[i].z = z;
		cloud_points[i].taken = false;
		
		center_x += x;
		center_y += y;
		center_z += z;
	}
	// Calculate mean point
	rootPoint[0] = center_x / n;
	rootPoint[1] = center_y / n;
	rootPoint[2] = center_z / n;

	// Read the scalars
	/*
	while (file >> str && str != "default");

	file >> str;
	printf("%.10lf\n",str.c_str());

	for (int i = 0; i < n; i++)
	{
		file >> scalar;
		cloud_points[i].scalar = scalar;
	}
	file.close();

	// Sort the points in ascending order of the scalar (activation time)
	qsort(cloud_points,n,sizeof(Point),compare);
	*/

	// DEBUG
	//printf("Root point = (%.10lf,%.10lf,%.10lf)\n",rootPoint[0],rootPoint[1],rootPoint[2]); 
}

void Miocardium::read_terminal_points (const std::string terminals_filename)
{
	printf("[Miocardium] Reading terminal points ...\n");	

	int n, m;
	double x, y, z;
	std::ifstream file;
	std::string str;

	file.open(terminals_filename.c_str());
	if (!file) 
	{
		printf("[-] ERROR! Cannot load terminals file '%s'!\n",terminals_filename.c_str());
		exit(EXIT_FAILURE);
	}

	file >> n >> m;
	this->num_terminal_points = n;
	this->terminal_points = new Point[n];

	Point *terminal_points = this->terminal_points;

	for (int i = 0; i < n; i++)
	{
		file >> x >> y >> z;

		terminal_points[i].x = x;
		terminal_points[i].y = y;
		terminal_points[i].z = z;
		terminal_points[i].taken = false;
	}

	file.close();	

	// DEBUG
	//this->the_miocardium->print();
	
}

// Discover the limits of the miocardium mesh on each coordinate (x,y,z)
void Miocardium::set_limits ()
{
	int n;
	double x, y, z;
	for (int i = 0; i < 3; i++)
	{	
		min_xyz[i] = 1.0e+10;
		max_xyz[i] = -1.0e+10;
	}

	n = this->num_cloud_points;
	// Pass through all the points of the mesh
	for (int i = 0; i < n; i++)
	{
		x = cloud_points[i].x; 
		y = cloud_points[i].y; 
		z = cloud_points[i].z;

		// Check the minimum values
		if (x < min_xyz[0]) min_xyz[0] = x;
		if (y < min_xyz[1]) min_xyz[1] = y;
		if (z < min_xyz[2]) min_xyz[2] = z;

		// Check the maximum values
		if (x > max_xyz[0]) max_xyz[0] = x;
		if (y > max_xyz[1]) max_xyz[1] = y;
		if (z > max_xyz[2]) max_xyz[2] = z;
	}
	printf("=============== LIMITS ================\n");
	printf("x - [%lf,%lf]\n",min_xyz[0],max_xyz[0]);
	printf("y - [%lf,%lf]\n",min_xyz[1],max_xyz[1]);
	printf("z - [%lf,%lf]\n",min_xyz[2],max_xyz[2]);
	printf("=======================================\n");
}


void Miocardium::print ()
{
	int n = this->num_cloud_points;
	Point *cloud_points = this->cloud_points;

	printf("-------------- Cloud points --------------\n");
	for (int i = 0; i < n; i++)
		printf("Point %d ---- %.10lf %.10lf %.10lf ---- AT = %.10lf\n",i,cloud_points[i].x,cloud_points[i].y,cloud_points[i].z,cloud_points[i].scalar);

	int m = this->num_terminal_points;
	Point *terminal_points = this->terminal_points;
	printf("-------------- Terminal points --------------\n");
	for (int i = 0; i < m; i++)
		printf("Terminal %d ---- %.10lf %.10lf %.10lf\n",i,terminal_points[i].x,terminal_points[i].y,terminal_points[i].z);
	
}

int compare (const void *a, const void *b)
{
	Point *pA = (Point*)a;
	Point *pB = (Point*)b;

	if (pA->scalar < pB->scalar)
      return -1;
   	else if (pA->scalar > pB->scalar)
      return 1;
   	else
      return 0;
}
