#include "../include/utils.h"

/*
void write_VTK_to_file (const double *sv, const double dx,\
                        const int Ncell, const int Nodes,\
			const int iter)
{
	FILE *file;
	int np, nl;
	char filename[50];

	np = Ncell;
	nl = Ncell-1;

	// Write the points
	sprintf(filename,"vtk/vm%d.vtk",iter);
	file = fopen(filename,"w+");
	fprintf(file,"# vtk DataFile Version 3.0\n");
	fprintf(file,"Monodomain FDM\n");
	fprintf(file,"ASCII\n");
	fprintf(file,"DATASET POLYDATA\n");
	fprintf(file,"POINTS %d float\n",np);
	for (int i = 0; i < Ncell; i++)
	{
		double x = dx*i;
		fprintf(file,"%g %g %g\n",x,0.0,0.0);
	}

	// Write the lines
	fprintf(file,"LINES %d %d\n",nl,nl*3);
	for (int i = 0; i < nl; i++)
	{
		fprintf(file,"2 %d %d\n",i,i+1);
	}

	// Write the transmembrane potential
	fprintf(file,"POINT_DATA %d\n",np);
	fprintf(file,"SCALARS vm float 1\n");
	fprintf(file,"LOOKUP_TABLE default\n");
	for (int i = 0; i < Ncell; i++)
	{
		fprintf(file,"%g\n",sv[i*Nodes]);
	}

	fclose(file);
}

void write_steady_state_to_file (const double *sv, const int Ncell, const int Nodes)
{
	FILE *file = fopen("output.sst","w+");
	for (int i = 0; i < Ncell; i++)
	{
		for (int j = 0; j < Nodes-1; j++)
		{
			fprintf(file,"%.10lf ",sv[i*Nodes+j]);
		}
		fprintf(file,"%.10lf\n",sv[i*Nodes+(Nodes-1)]);
	}
	fclose(file);
}

void write_plot_data (std::ofstream files[], const double t, const double *sv,\
			const int Ncell, const int Nodes,\
			const int ids[])
{
	for (int i = 0; i < 5; i++)
	{
		int id = ids[i];
		files[i] << std::setprecision(10) << t << " " << sv[id*Nodes] << " " << sv[id*Nodes+1] << std::endl;
		//fprintf(files[i],"%.10lf %.10lf\n",t,sv[id*Nodes],sv[id*Nodes+1]);
	}
}

void print_stimulus (const double *stim_current, const int Ncell, const double dx)
{
	for (int i = 0; i < Ncell; i++)
	{
		double x = i*dx;
		printf("Cell at x = %.10lf ---- Istim = %g\n",x,stim_current[i]);
	}
}

void print_state_vector (const double *sv, const int Ncell, const int Nodes)
{
	printf("%s\n",PRINT_LINE);
	for (int i = 0; i < Ncell; i++)
	{
		printf("Cell %d --> ",i);
		for (int j = 0; j < Nodes-1; j++)
		{
			printf("%.10lf ",sv[i*Nodes+j]);
		}
		printf("%.10lf\n",sv[i*Nodes+(Nodes-1)]);
	}
	printf("%s\n",PRINT_LINE);
}

void print_configuration_parameters(const double dx, const double dt, const double tmax, const double lmax,\
					const int Ncell, const int Niter, const int Nodes,\
					const int plot_cell_ids[], const int print_rate, const int sst_rate)
{
	int i;

	printf("%s\n",PRINT_LINE);
	printf("[Utils] dx = %.10lf cm\n",dx);
	printf("[Utils] dt = %.10lf ms\n",dt);
	printf("[Utils] tmax = %.10lf ms\n",tmax);
	printf("[Utils] lmax = %.10lf cm\n",lmax);
	printf("[Utils] Number of cells = %d\n",Ncell);
	printf("[Utils] Number of iterations = %d\n",Niter);
	printf("[Utils] Number of ODE equations = %d\n",Nodes);
	printf("[Utils] Print rate = %d\n",print_rate);
	printf("[Utils] Sst rate = %d\n",sst_rate);
	printf("[Utils] Plot cell ids = ");
	for (i = 0; i < 4; i++)
		printf("%d ",plot_cell_ids[i]);
	printf("%d\n",plot_cell_ids[i]);
	printf("%s\n",PRINT_LINE);
}


void print_progress (int iter, int max_iter)
{
    double progress = iter / (double)max_iter;
    
    std::cout << "Progress: " << int(progress * 100.0) << " %\r";
    std::cout.flush();
}

void configure_plot_cells (std::ofstream files[], int plot_cell_ids[], const double lmax, const double dx)
{
	int offset = nearbyint(lmax / dx) / 5;
	for (int i = 0; i < 5; i++)
	{
		std::stringstream filename;
		plot_cell_ids[i] = i*offset;	
	
		filename << "output/sv-" << plot_cell_ids[i] << ".dat";
		files[i].open(filename.str().c_str());		
	}
}

void close_plot_files (std::ofstream files[])
{
	for (int i = 0; i < 5; i++)
		files[i].close();
}
*/

void Usage (const char pname[])
{
	printf("%s\n",PRINT_LINE);
	printf("Usage:> %s <input_config_file>\n",pname);
	printf("\t<input_config_file> = Input file with the parameters values for the simulation\n");
	printf("%s\n",PRINT_LINE);	
}
