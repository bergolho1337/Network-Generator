// Author: Lucas Berg
// =============================================================================================
// Program that reads 3 files and merge into a single '.vtk' file.
// 
// - Input File 1 -> Bundle-His.pkje
// - Input File 2 -> LV-Purkinje.pkje
// - Input File 3 -> RV-Purkinje.pkje
//
// Pre-Requisites:
//  - GCC/G++ 
//  - CMake
//  - Valgrind (for memory leak check ...)
//
//  Features of this program:
//
//   1) Read 3 Purkinje networks in '.vkt' format and merge them together into another '.vtk' file
//   2) Write the merged network in '.vtk' format  
//
// How to compile:
//  $ ./recompile_project.sh
// =============================================================================================

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdint>
#include <cmath>
#include <vector>

uint32_t total_number_nodes;
uint32_t total_number_edges;

struct node
{
    uint32_t id;
    double x, y, z;
    std::vector< std::pair<uint32_t,double> > edges;
};

void new_purkinje_network(const char filename[], std::vector< struct node > &pk)
{
    FILE *file = fopen(filename,"r");

    char str[200];
    while (fscanf(file,"%s",str) != EOF)
    {
        if (strcmp(str,"POINTS") == 0) break;
    }
    int num_nodes;
    fscanf(file,"%d %s",&num_nodes,str);

    for (int i = 0; i < num_nodes; i++)
    {
        double x, y, z;
        fscanf(file,"%lf %lf %lf",&x,&y,&z);

        struct node n;
        n.id = i;
        n.x = x;
        n.y = y;
        n.z = z;
        pk.push_back(n);
    }

    while (fscanf(file,"%s",str) != EOF)
    {
        if (strcmp(str,"LINES") == 0) break;
    }    
    int num_edges, trash;
    std::vector< std::pair<int,int> > lines;
    fscanf(file,"%d %d",&num_edges,&trash);
    for (int i = 0; i < num_edges; i++)
    {
        int src, dest;
        fscanf(file,"%d %d %d",&trash,&src,&dest);

        pk[src].edges.push_back( std::make_pair(dest,0) );
        lines.push_back( std::make_pair(src,dest) );
    }
    while (fscanf(file,"%s",str) != EOF)
    {
        if (strcmp(str,"default") == 0) break;
    }
    
    for (int i = 0; i < num_edges; i++)
    {
        double value;
        fscanf(file,"%lf",&value);

        int src = lines[i].first;
        int dest = lines[i].second;
        for (uint32_t j = 0; j < pk[src].edges.size(); j++)
        {
            if (pk[src].edges[j].first == dest)
                pk[src].edges[j].second = value;
        }
    }

    fclose(file);
} 

void print_network (std::vector< struct node > pk)
{
    for (uint32_t i = 0; i < pk.size(); i++)
    {
        printf("|| %d %g %g %g || ",pk[i].id,pk[i].x,pk[i].y,pk[i].z);
        for (uint32_t j = 0; j < pk[i].edges.size(); j++)
        {
            printf(" --> || %d ||",pk[i].edges[j].first);
        }
        printf("\n");
    }
}

uint32_t get_root_index (std::vector< struct node > lv_network, const double root[3])
{
    double min_dist = __DBL_MAX__;
    uint32_t min_index = 0;

    for (uint32_t i = 0; i < lv_network.size(); i++)
    {
        double p1[3];
        p1[0] = lv_network[i].x;
        p1[1] = lv_network[i].y;
        p1[2] = lv_network[i].z;

        double dist = sqrt(pow(p1[0]-root[0],2) + pow(p1[1]-root[1],2) + pow(p1[2]-root[2],2));
        if (dist < min_dist)
        {
            min_dist = dist;
            min_index = i;
        }
    }
    return min_index;
}

void merge_purkinje_networks (std::vector< struct node > bundle_his, std::vector< struct node > lv_network, std::vector< struct node > rv_network,\
                            std::vector< struct node > &new_network)
{
    uint32_t offset = 0;
    for (uint32_t i = 0; i < bundle_his.size(); i++)
    {
        struct node n;
        n.id = i;
        n.x = bundle_his[i].x;
        n.y = bundle_his[i].y;
        n.z = bundle_his[i].z;

        new_network.push_back(n);
    }
    offset = bundle_his.size();
    uint32_t first_lv_index = offset;
    for (uint32_t i = 0; i < lv_network.size(); i++)
    {
        struct node n;
        n.id = i + offset;
        n.x = lv_network[i].x;
        n.y = lv_network[i].y;
        n.z = lv_network[i].z;

        new_network.push_back(n);
    }
    offset += lv_network.size();
    uint32_t first_rv_index = offset;
    for (uint32_t i = 0; i < rv_network.size(); i++)
    {
        struct node n;
        n.id = i + offset;
        n.x = rv_network[i].x;
        n.y = rv_network[i].y;
        n.z = rv_network[i].z;

        new_network.push_back(n);
    }
    total_number_nodes = bundle_his.size() + lv_network.size() + rv_network.size();

    // Set edges for the Bundle of His
    offset = 0;
    total_number_edges = 0;
    for (uint32_t i = 0; i < bundle_his.size(); i++)
    {
        int src = i;
        for (uint32_t j = 0; j < bundle_his[i].edges.size(); j++)
        {
            uint32_t dest = bundle_his[i].edges[j].first;
            double radius = bundle_his[i].edges[j].second;
            new_network[src].edges.push_back( std::make_pair(dest,radius) );
            total_number_edges++;
        }
    }

    // Set edges for the LV
    offset = bundle_his.size();
    for (uint32_t i = 0; i < lv_network.size(); i++)
    {
        uint32_t src = i;
        for (uint32_t j = 0; j < lv_network[i].edges.size(); j++)
        {
            uint32_t dest = lv_network[i].edges[j].first;
            double radius = lv_network[i].edges[j].second;
            new_network[src+offset].edges.push_back( std::make_pair(dest+offset,radius) );
            total_number_edges++;
        }
    }
    // Set edges for the RV
    offset += lv_network.size();
    for (uint32_t i = 0; i < rv_network.size(); i++)
    {
        uint32_t src = i;
        for (uint32_t j = 0; j < rv_network[i].edges.size(); j++)
        {
            uint32_t dest = rv_network[i].edges[j].first;
            double radius = rv_network[i].edges[j].second;
            new_network[src+offset].edges.push_back( std::make_pair(dest+offset,radius) );
            total_number_edges++;
        }
    }

    new_network[1].edges.push_back( std::make_pair(first_lv_index,100.0) );
    new_network[1].edges.push_back( std::make_pair(first_rv_index,100.0) );
    total_number_edges += 2;
}

void write_network_to_vtk(const char filename[], std::vector< struct node > new_network)
{
    FILE *file = fopen(filename,"w+");

    fprintf(file,"# vtk DataFile Version 3.0\n");
    fprintf(file,"Tree\n");
    fprintf(file,"ASCII\n");
    fprintf(file,"DATASET POLYDATA\n");
    fprintf(file,"POINTS %lu float\n",new_network.size());
    for (uint32_t i = 0; i < new_network.size(); i++)
        fprintf(file,"%g %g %g\n",new_network[i].x,new_network[i].y,new_network[i].z);
    fprintf(file,"LINES %u %u\n",total_number_edges,total_number_edges*3);
    for (uint32_t i = 0; i < new_network.size(); i++)
    {
        for (uint32_t j = 0; j < new_network[i].edges.size(); j++)
        {
            fprintf(file,"2 %u %u\n",i,new_network[i].edges[j].first);
        }
    }
    
    fprintf(file,"CELL_DATA %u\n",total_number_edges);
    fprintf(file,"SCALARS radius float\n");
    fprintf(file,"LOOKUP_TABLE default\n");
    for (uint32_t i = 0; i < new_network.size(); i++)
    {
        for (uint32_t j = 0; j < new_network[i].edges.size(); j++)
        {
            fprintf(file,"%g\n",new_network[i].edges[j].second);
        }
    }

    //fprintf(file,"VERTICES %u %u\n",new_network.size(),new_network.size()*2);
    //for (uint32_t i = 0; i < new_network.size(); i++)
    //    fprintf(file,"1 %u\n",i);


    fclose(file);
}

int main (int argc, char *argv[])
{
    if (argc-1 != 4)
    {
        printf("Usage:> %s <input_filename_1> <input_filename_2> <input_filename_3> <output_filename>\n",argv[0]);
        exit(EXIT_FAILURE);
    }

    char *bundle_his_filename = argv[1];
    char *lv_network_filename = argv[2];
    char *rv_network_filename = argv[3];
    char *lvrv_network_filename = argv[4];

    std::vector< struct node > bundle_his;
    new_purkinje_network(bundle_his_filename,bundle_his);
    //print_network(bundle_his);

    std::vector< struct node > lv_network;
    new_purkinje_network(lv_network_filename,lv_network);
    //print_network(lv_network);

    std::vector< struct node > rv_network;
    new_purkinje_network(rv_network_filename,rv_network);
    //print_network(rv_network);

    std::vector< struct node > new_network;
    merge_purkinje_networks(bundle_his,lv_network,rv_network,new_network);
    //print_network(new_network);

    write_network_to_vtk(lvrv_network_filename,new_network);

    
    return 0;
}