// Author: Lucas Berg
// Program that receive a Purkinje network written on .vtk format and returns the duplicate points on it.

#include <iostream>
#include <vector>
#include <algorithm>
#include <map>

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>

using namespace std;

// ====================================================================================================================
// CONSTANTS AND MACROS
#define PRINT_LINE "--------------------------------------------------------------------------------------------"
// ====================================================================================================================

// ====================================================================================================================
// CLASSES
class Point_3D
{
public:
    double x, y, z;
public:
    Point_3D (const double x, const double y, const double z)
    {
        this->x = x;
        this->y = y;
        this->z = z;
    }

    bool operator <(const Point_3D& pt) const
    {
        return (this->x < pt.x) || \
               ((!(pt.x < this->x)) && (this->y < pt.y)) || \
               ((!(pt.x < this->x)) && (!(pt.y < this->y)) && (this->z < pt.z));
    }
    friend std::ostream& operator<<(std::ostream& os, const Point_3D& obj);
};

class Line
{
public:
    uint32_t src, dest;
public:
    Line (const uint32_t src, const uint32_t dest)
    {
        this->src = src;
        this->dest = dest;
    }
};
// ====================================================================================================================

int main (int argc, char *argv[])
{
  if (argc-1 != 1)
  {
      printf("%s\n",PRINT_LINE);
      printf("Usage:> %s <input_filename>\n",argv[0]);
      printf("%s\n",PRINT_LINE);
      printf("<input_filename> = Input filename on .vtk format\n");
      printf("%s\n",PRINT_LINE);
      exit(EXIT_FAILURE);
  }

  char *filename = argv[1];

  std::vector<Point_3D> points;
  std::vector<Line> lines;

  uint32_t unique_num_points = 0;
  std::map<Point_3D,uint32_t> points_map;
  std::map<uint32_t,uint32_t> duplicates_map;
  std::map<Point_3D,uint32_t>::iterator it;
  std::map<uint32_t,uint32_t>::iterator it2;

  // ==============================================================================================================
  // READING
  FILE *file = fopen(filename,"r");
  if (!file)
  {
      fprintf(stderr,"[-] ERROR! Cannot open file '%s'\n",filename);
      exit(EXIT_FAILURE);
  }

  uint32_t total_nodes, total_edges;
  char str[100], trash[100];
  double pos[3];
  uint32_t edge[3];

  // Read the nodes
  while (fscanf(file,"%s",str) != EOF)
      if (strcmp(str,"POINTS") == 0) break;
  fscanf(file,"%u",&total_nodes);
  fscanf(file,"%s",trash);

  for (uint32_t i = 0; i < total_nodes; i++)
  {
      fscanf(file,"%lf %lf %lf",&pos[0],&pos[1],&pos[2]);

      Point_3D point(pos[0],pos[1],pos[2]);

      points.push_back(point);

      it = points_map.find(point);

      if (it == points_map.end())
      {
          points_map.insert(std::pair<Point_3D,uint32_t>(point,unique_num_points));

          unique_num_points++;
      }
      else
      {
        uint32_t original = it->second;
        uint32_t duplicate = i;

        duplicates_map.insert(std::pair<uint32_t,uint32_t>(duplicate,original));
        printf("Current point %u -- [!] Duplicate point -- (%lf %lf %lf) --> (%u)\n",duplicate,pos[0],pos[1],pos[2],original);
      }
  }

  // Read the edges
  while (fscanf(file,"%s",str) != EOF)
      if (strcmp(str,"LINES") == 0) break;
  fscanf(file,"%u",&total_edges);
  fscanf(file,"%s",trash);
  for (uint32_t i = 0; i < total_edges; i++)
  {
      fscanf(file,"%u %u %u",&edge[0],&edge[1],&edge[2]);

      Line line(edge[1],edge[2]);

      lines.push_back(line);
  }

  fclose(file);

/*
  // ==============================================================================================================
  // ELIMINATE DUPLICATES
  for (it2 = duplicates_map.begin(); it2 != duplicates_map.end(); ++it2)
  {
    uint32_t duplicate_index = it2->first;
    uint32_t original_index = it2->second;

    // Eliminate every line that contains a duplicate points either on 'src' or 'dest'
    for (uint32_t i = 0; i < lines.size(); i++)
    {
      uint32_t src = lines[i].src;
      uint32_t dest = lines[i].dest;

      if (src == duplicate_index || dest == duplicate_index)
      {
        lines.erase(lines.begin()+i);
        i--;
      }
    }
  }

  // Eliminate duplicate point
  uint32_t prev_duplicate_index = points.size()+1;
  for (it2 = duplicates_map.begin(); it2 != duplicates_map.end(); ++it2)
  {
    uint32_t duplicate_index = it2->first;
    uint32_t original_index = it2->second;

    if (duplicate_index > prev_duplicate_index)
      duplicate_index--;

    points.erase(points.begin()+duplicate_index);

    for (uint32_t i = 0; i < lines.size(); i++)
    {

      if (lines[i].src > duplicate_index)
      {
        lines[i].src--;
      }
      if (lines[i].dest > duplicate_index)
      {
        lines[i].dest--;
      }
    }

    prev_duplicate_index = duplicate_index;
  }


  // ==============================================================================================================
  // WRITING
  file = fopen("outputs/processed.vtk","w+");

  fprintf(file,"# vtk DataFile Version 3.0\n");
  fprintf(file,"Purkinje\n");
  fprintf(file,"ASCII\n");
  fprintf(file,"DATASET POLYDATA\n");
  fprintf(file,"POINTS %lu float\n",points.size());
  for (uint32_t i = 0; i < points.size(); i++)
      fprintf(file,"%g %g %g\n",points[i].x,points[i].y,points[i].z);
  fprintf(file,"LINES %lu %lu\n",lines.size(),lines.size()*3);
  for (uint32_t i = 0; i < lines.size(); i++)
      fprintf(file,"2 %u %u\n",lines[i].src,lines[i].dest);

  fclose(file);
*/
  return 0;
}
