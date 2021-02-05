// ********************************************************************************************************************
// 	POINT REMAPPER - GUIDED EDITION
//	Author: Lucas Berg
//
//	Input arguments:
//		- A new version of the Point-Remapper that remap the points from a cloud of points using a given
//	Purkinje network as reference.
//	Outputs:
//		- A PTS file with the remapped points following the input Purkinje network BFS search path.
//
// ********************************************************************************************************************

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <fstream>

#include "reader/reader.h"

using namespace std;

void usage (const char pname[])
{
    printf("%s\n",PRINT_LINE);
    printf("Usage:> %s <network_file> <cloud_point_file>\n",pname);
    printf("\t<network_file> = Input file with the graph configuration\n");
    printf("\t<cloud_point_file> = Input file with the cloud points\n");
    printf("%s\n",PRINT_LINE);
}


int main (int argc, char *argv[])
{
	if (argc-1 != 2)
	{
		usage(argv[0]);
		exit(EXIT_FAILURE);
	}

	Reader *reader = new Reader(argc,argv);

	reader->remap_points_using_graph();

	delete reader;

	return 0;  
}
