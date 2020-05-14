// Author: Lucas Berg
// Cable equation from Effects of Diameter ....
// v = sqrt( (G_i * d) / (4 * C_f * tau_f) )

#include <iostream>
#include <vector>
#include <algorithm>

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>

using namespace std;

// ====================================================================================================================
// CONSTANTS AND MACROS
#define PRINT_LINE "--------------------------------------------------------------------------------------------"
#define UM_TO_CM 0.0001
#define UM_TO_M 0.000001
#define MS_TO_S 0.001
#define S_TO_MS 1000.0
#define CM_S_TO_M_S 0.01

const uint32_t NPOINTS = 200;
const double G_I = 7.9;                     // mmho/cm
const double C_F = 3.4;                     // uF/cm^2
const double TAU_F = 0.1;                   // ms

const double MIN_DIAMETER = 100.0;             // um
const double MAX_DIAMETER = 200.0;             // um
// ====================================================================================================================

// m/s
double calc_velocity (const double diameter)
{
  return powf( (G_I*diameter)/(4.0*C_F*TAU_F) , 0.5 ) * 0.1;
}

// s
double calc_activation_time (const double diameter, const double length)
{
  return (length*UM_TO_M) / calc_velocity(diameter);
}

/*
void write_datafile (const double Rm, const double Rc, const double Cm, const double c)
{
  FILE *file_s = fopen("output/propagation_velocity.txt","w+");
  FILE *file_at = fopen("output/activation_time.txt","w+");

  double delta = (MAX_RADIUS - MIN_RADIUS) / (double)NPOINTS;

  for (uint32_t i = 0; i < NPOINTS; i++)
  {
    double radius = MIN_RADIUS + delta*i;
    double diameter = 2.0*radius;

    double lambda_m = calc_lambda_m(Rm,Rc,diameter);
    double tau_m = calc_tau_m(Rm,Cm);
    double s = calc_velocity(lambda_m,tau_m,c)*10.0;    // mm/ms
    double at = CABLE_LENGTH / s;                       // ms

    fprintf(file_s,"%g %g\n",diameter,s);
    fprintf(file_at,"%g %g\n",diameter,at);
  }
  fclose(file_s);
  fclose(file_at);
}
*/

class Node
{
public:
  double x, y, z;
  vector<uint32_t> edges;
public:
  Node () {}
};

void read_graph_from_vtk (const char filename[], vector<Node> &nodes)
{
  FILE *file = fopen(filename,"r");

  char str[200];
  while (fscanf(file,"%s",str) != EOF)
  {
    if (strcmp(str,"POINTS") == 0)
      break;
  }

  uint32_t num_nodes;
  fscanf(file,"%u",&num_nodes);
  fscanf(file,"%s",str);
  for (uint32_t i = 0; i < num_nodes; i++)
  {
    double pos[3];
    fscanf(file,"%lf %lf %lf",&pos[0],&pos[1],&pos[2]);

    Node node;
    node.x = pos[0];
    node.y = pos[1];
    node.z = pos[2];
    nodes.push_back(node);
  }

  while (fscanf(file,"%s",str) != EOF)
  {
    if (strcmp(str,"LINES") == 0)
      break;
  }

  uint32_t num_edges;
  fscanf(file,"%u",&num_edges);
  fscanf(file,"%s",str);
  uint32_t src, dest, trash;
  for (uint32_t i = 0; i < num_edges; i++)
  {
    fscanf(file,"%u %u %u",&trash,&src,&dest);
    nodes[src].edges.push_back(dest);
  }

  fclose(file);
}

void print_graph (vector<Node> nodes)
{
  for (uint32_t i = 0; i < nodes.size(); i++)
  {
    printf("|| %u %g %g %g || --> ",i,nodes[i].x,nodes[i].y,nodes[i].z);
    for (uint32_t j = 0; j < nodes[i].edges.size(); j++)
      printf("|| %u || --> ",nodes[i].edges[j]);
    printf("\n");
  }
}

void print_propagation_velocity ()
{
  double delta = (MAX_DIAMETER - MIN_DIAMETER) / (double)NPOINTS;

  for (uint32_t i = 0; i < NPOINTS; i++)
  {
    double diameter = MIN_DIAMETER + delta*i;

    double v = calc_velocity(diameter);                       // m/s
    printf("d = %g -- v = %g\n",diameter,v);
  }
}

int main (int argc, char *argv[])
{
  vector<Node> nodes;
  read_graph_from_vtk("inputs/elizabeth_largest_segment_LV.vtk",nodes);
  //print_graph(nodes);

  print_propagation_velocity();

  double total_activation_time = 0.0;
  double total_length = 0.0;
  for (uint32_t i = 0; i < nodes.size(); i++)
  {
    double pos1[3];
    pos1[0] = nodes[i].x;
    pos1[1] = nodes[i].y;
    pos1[2] = nodes[i].z;

    for (uint32_t j = 0; j < nodes[i].edges.size(); j++)
    {
      double pos2[3];
      uint32_t dest_index = nodes[i].edges[j];

      pos2[0] = nodes[dest_index].x;
      pos2[1] = nodes[dest_index].y;
      pos2[2] = nodes[dest_index].z;

      double length = sqrt( powf(pos2[0]-pos1[0],2) + powf(pos2[1]-pos1[1],2) + powf(pos2[2]-pos1[2],2) );

      double at = calc_activation_time(100.0,length);
      total_length += length;
      total_activation_time += at;
    }
  }
  printf("Activation time = %g -- Length = %g -- Velocity = %g\n",total_activation_time*S_TO_MS,total_length*UM_TO_CM,calc_velocity(100));

  return 0;
}
