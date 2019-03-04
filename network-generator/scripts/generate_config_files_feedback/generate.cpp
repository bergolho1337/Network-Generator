#include <iostream>
#include <cstdio>

using namespace std;

const int MAX_SIZE = 500;

void write_sst_config_file (const int period, const int start_period, const int step_period)
{
    char filename[MAX_SIZE];
    sprintf(filename,"files/sst_cable_5cm_%dms.ini",period);
    FILE *file = fopen(filename,"w+");

    fprintf(file,"[main]\n");
    fprintf(file,"dx = 0.01\n");
    fprintf(file,"dt = 0.01\n");
    fprintf(file,"tmax = %d\n",20*period);
    fprintf(file,"lmax = 5.0\n");
    fprintf(file,"print_rate = 100\n");
    fprintf(file,"sst_rate = %d\n",2000*period);
    fprintf(file,"use_steady_state = no\n");
    fprintf(file,"sst_filename = teste.sst\n");
    /*
    if (period == start_period)
    {
	fprintf(file,"use_steady_state = no\n");
        fprintf(file,"sst_filename = teste.sst\n");
    }
    else
    {
	fprintf(file,"use_steady_state = yes\n");
        fprintf(file,"sst_filename = steady_state/cable-5cm-%dms.sst\n",period+step_period);
    }
    */
    fprintf(file,"\n\n");
    fprintf(file,"[stimulus]\n");
    fprintf(file,"stim_start = 0.0\n");
    fprintf(file,"stim_duration = 5.0\n");
    fprintf(file,"start_period = %d\n",period);
    fprintf(file,"end_period = %d\n",period);
    fprintf(file,"period_step = 100\n");
    fprintf(file,"n_cycles = 20\n");
    
    fclose(file);
}

void write_simulation_config_file (const int period, const int start_period)
{

    char filename[MAX_SIZE];
    sprintf(filename,"files/simple_cable_5cm_%dms.ini",period);
    FILE *file = fopen(filename,"w+");

    fprintf(file,"[main]\n");
    fprintf(file,"dx = 0.01\n");
    fprintf(file,"dt = 0.01\n");
    fprintf(file,"tmax = %d\n",3*period);
    fprintf(file,"lmax = 5.0\n");
    fprintf(file,"print_rate = 100\n");
    fprintf(file,"sst_rate = %d\n",2000*period);
    fprintf(file,"use_steady_state = yes\n");
    fprintf(file,"sst_filename = steady_state/cable-5cm-%dms.sst\n",period);
    fprintf(file,"\n\n");
    fprintf(file,"[stimulus]\n");
    fprintf(file,"stim_start = 0.0\n");
    fprintf(file,"stim_duration = 5.0\n");
    fprintf(file,"start_period = %d\n",period);
    fprintf(file,"end_period = %d\n",period);
    fprintf(file,"period_step = 100\n");
    fprintf(file,"n_cycles = 20\n");
    
    fclose(file);
}

void write_config_files (const int period, const int start_period, const int step_period)
{
    write_sst_config_file(period,start_period,step_period);
    write_simulation_config_file(period,start_period);
}

int main ()
{
    int start_period = 300;
    int end_period = 100;
    int period_step = 5;
    for (int period = start_period; period >= end_period; period -= period_step)
    {
        printf("Writing %dms files ...\n",period);
        write_config_files(period,start_period,period_step);
    }
    return 0;
}
