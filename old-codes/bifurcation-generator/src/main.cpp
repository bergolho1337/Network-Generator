// Author: Lucas Berg

#include <iostream>
#include <vector>
#include <algorithm>

#include <cstdio>
#include <cstdlib>
#include <cmath>

using namespace std;

// ====================================================================================================================
// CONSTANTS AND MACROS
#define PRINT_LINE "--------------------------------------------------------------------------------------------"
// ====================================================================================================================

struct point3d
{
  uint32_t id;
  double x, y, z;
};

struct segment
{
  uint32_t id;
  struct point3d *src;
  struct point3d *dest;
};

struct network
{
  uint32_t number_of_points;
  uint32_t number_of_segments;
  struct point3d *the_points;
  struct segment *the_segments;
};

struct network* new_network ()
{
    struct network *result = (struct network*)malloc(sizeof(struct network));

    result->number_of_points = 0;
    result->number_of_segments = 0;
    result->the_points = NULL;
    result->the_segments = NULL;

    return result;
}

void make_root (struct network *the_network, const double L)
{
  // Set the 2 first points and the main segment of the network
  the_network->the_points[0].id = 0;
  the_network->the_points[0].x = 0;
  the_network->the_points[0].y = 0;
  the_network->the_points[0].z = 0;

  the_network->the_points[1].id = 1;
  the_network->the_points[1].x = L;
  the_network->the_points[1].y = 0;
  the_network->the_points[1].z = 0;

  the_network->the_segments[0].id = 0;
  the_network->the_segments[0].src = &the_network->the_points[0];
  the_network->the_segments[0].dest = &the_network->the_points[1];
}

void make_bifurcation (struct network *the_network, const double L, const uint32_t num_biff)
{
  struct point3d *the_points = the_network->the_points;
  struct segment *the_segments = the_network->the_segments;

  // Get a reference to coordinate of teh bifurcation point
  double x = the_points[1].x;
  double y = the_points[1].y;
  double z = the_points[1].z;

  double angle = M_PI / 4.0;
  double half_angle = angle / 2.0;

  switch (num_biff)
  {
    case 1: {
              the_points[2].id = 2;
              the_points[2].x = 2*L;
              the_points[2].y = 0;
              the_points[2].z = 0;

              the_segments[1].id = 1;
              the_segments[1].src = &the_points[1];
              the_segments[1].dest = &the_points[2];

              break;
            }
    case 2: {
              the_points[2].id = 2;
              the_points[2].x = x + cos(angle)*L;
              the_points[2].y = y + sin(angle)*L;
              the_points[2].z = 0;

              the_segments[1].id = 1;
              the_segments[1].src = &the_points[1];
              the_segments[1].dest = &the_points[2];

              the_points[3].id = 3;
              the_points[3].x = x + cos(-angle)*L;
              the_points[3].y = y + sin(-angle)*L;
              the_points[3].z = 0;

              the_segments[2].id = 2;
              the_segments[2].src = &the_points[1];
              the_segments[2].dest = &the_points[3];

              break;
            }
    case 3: {
              the_points[2].id = 2;
              the_points[2].x = x + cos(angle)*L;
              the_points[2].y = y + sin(angle)*L;
              the_points[2].z = 0;

              the_segments[1].id = 1;
              the_segments[1].src = &the_points[1];
              the_segments[1].dest = &the_points[2];

              the_points[3].id = 3;
              the_points[3].x = x + cos(-angle)*L;
              the_points[3].y = y + sin(-angle)*L;
              the_points[3].z = 0;

              the_segments[2].id = 2;
              the_segments[2].src = &the_points[1];
              the_segments[2].dest = &the_points[3];

              the_points[4].id = 4;
              the_points[4].x = x + L;
              the_points[4].y = 0;
              the_points[4].z = 0;

              the_segments[3].id = 3;
              the_segments[3].src = &the_points[1];
              the_segments[3].dest = &the_points[4];

              break;
            }
      case 4: {
                the_points[2].id = 2;
                the_points[2].x = x + cos(angle)*L;
                the_points[2].y = y + sin(angle)*L;
                the_points[2].z = 0;

                the_segments[1].id = 1;
                the_segments[1].src = &the_points[1];
                the_segments[1].dest = &the_points[2];

                the_points[3].id = 3;
                the_points[3].x = x + cos(half_angle)*L;
                the_points[3].y = y + sin(half_angle)*L;
                the_points[3].z = 0;

                the_segments[2].id = 2;
                the_segments[2].src = &the_points[1];
                the_segments[2].dest = &the_points[3];

                the_points[4].id = 4;
                the_points[4].x = x + cos(-half_angle)*L;
                the_points[4].y = y + sin(-half_angle)*L;
                the_points[4].z = 0;

                the_segments[3].id = 3;
                the_segments[3].src = &the_points[1];
                the_segments[3].dest = &the_points[4];

                the_points[5].id = 5;
                the_points[5].x = x + cos(-angle)*L;
                the_points[5].y = y + sin(-angle)*L;
                the_points[5].z = 0;

                the_segments[4].id = 4;
                the_segments[4].src = &the_points[1];
                the_segments[4].dest = &the_points[5];

                break;
              }
  }
}

void build_network (struct network *the_network, const double L, const uint32_t num_biff)
{
  // Initialize and allocate memory
  the_network->number_of_points = num_biff + 2;
  the_network->number_of_segments = the_network->number_of_points - 1;
  the_network->the_points = (struct point3d*)malloc(sizeof(struct point3d)*the_network->number_of_points);
  the_network->the_segments = (struct segment*)malloc(sizeof(struct segment)*the_network->number_of_segments);

  make_root(the_network,L);

  make_bifurcation(the_network,L,num_biff);
}

void free_network (struct network *the_network)
{
  if (the_network->the_points)
    free(the_network->the_points);

  if (the_network->the_segments)
    free(the_network->the_segments);

  free(the_network);
}

void print_network (struct network *the_network)
{
  uint32_t number_of_points = the_network->number_of_points;
  uint32_t number_of_segments = the_network->number_of_segments;

  for (uint32_t i = 0; i < number_of_points; i++)
  {
    printf("Point %u = (%g,%g,%g)\n",the_network->the_points[i].id,the_network->the_points[i].x,the_network->the_points[i].y,the_network->the_points[i].z);
  }

  for (uint32_t i = 0; i < number_of_segments; i++)
  {
    printf("Segment %u = (%u,%u)\n",the_network->the_segments[i].id,the_network->the_segments[i].src->id,the_network->the_segments[i].dest->id);
  }
}

void write_network (struct network *the_network)
{
  FILE *file = fopen("output/network.vtk","w+");

  fprintf(file,"# vtk DataFile Version 3.0\n");
  fprintf(file,"Bifurcation network\n");
  fprintf(file,"ASCII\n");
  fprintf(file,"DATASET POLYDATA\n");
  fprintf(file,"POINTS %u float\n",the_network->number_of_points);
  for (uint32_t i = 0; i < the_network->number_of_points; i++)
    fprintf(file,"%g %g %g\n",the_network->the_points[i].x,the_network->the_points[i].y,the_network->the_points[i].z);
  fprintf(file,"LINES %u %u\n",the_network->number_of_segments,the_network->number_of_segments*3);
  for (uint32_t i = 0; i < the_network->number_of_segments; i++)
    fprintf(file,"2 %u %u\n",the_network->the_segments[i].src->id,the_network->the_segments[i].dest->id);

  fclose(file);
}

int main (int argc, char *argv[])
{
    if (argc-1 != 2)
    {
        printf("%s\n",PRINT_LINE);
        printf("Usage:> %s <L> <number_of_bifurcations>\n",argv[0]);
        printf("%s\n",PRINT_LINE);
        printf("<L> = Size of each segment {cm}\n");
        printf("<number_of_bifurcations> = Number of bifurcations\n");
        printf("%s\n",PRINT_LINE);
        exit(EXIT_FAILURE);
    }

    double L = atof(argv[1]);
    uint32_t number_of_bifurcations = atoi(argv[2]);

    struct network *the_network = new_network();

    build_network(the_network,L,number_of_bifurcations);
    //print_network(the_network);

    write_network(the_network);

    free_network(the_network);

    return 0;
}
