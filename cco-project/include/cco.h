#ifndef CCO_H
#define CCO_H

#include <cmath>
#include <ctime>
#include <vector>
#include <algorithm>

#include <vtkRegularPolygonSource.h>
#include <vtkXMLPolyDataWriter.h>


#include "options.h"

// Constants and Macros
const int ROOT_STICKER = 1;
const int TERMINAL_STICKER = 2;
const int BIFURCATION_STICKER = 3;
const int NIL = -1;
const double ETA = 3.6e-03;         // Viscosity of blood (Pa.s)
const int N_toss = 200;             // Number of tosses for a new terminal

#define PRINT_LINE "------------------------------------------------------------------------------"

using namespace std;

class Point
{
public:
    unsigned int id;
    double x, y, z;

public:
    Point ();
    Point (const int id, const double x, const double y, const double z);
};

class Segment
{
public:
    int type;

    double beta_l;
    double beta_r;
    double Q;
    double p;
    double radius;
    double resistance;          // Relative resistance
    double length;
    int ndist;

    unsigned int src;           // Source Point index
    unsigned int dest;          // Destination Point index
    //Point *p1;
    //Point *p2;

    int left;          // Left offspring Segment index
    int right;         // Right offspring Segment index
    int parent;        // Parent offspring Segment index
    //Segment *left;
    //Segment *right;
    //Segment *parent;
public:
    Segment ();
    Segment (Point *p1, Point *p2,\
            int left, int right, int parent,\
            const double Q, const double p);

    void add_offspring (Segment *new_segment);
    double calc_dproj (const Point p);
    double calc_dortho (const Point p);
    double calc_dend (const Point p);
private:
    int get_segment_type ();
};

class CCO_Network
{
private:
    int num_terminals;

    int N_term;
    double Q_perf;
    double p_perf;
    double r_perf;
    double r_supp;
    
    vector<Point> points;
    vector<Segment> segments;
    Point center;

public:
    CCO_Network ();
    CCO_Network (User_Options *options);
    void grow_tree ();
    void make_root ();
    void test1 ();
    void test2 ();
    void test3 ();
    void test4 ();

    void generate_new_terminal ();
    void build_segment (double new_pos[]);
    void build_segment (const unsigned int j);
    void build_segment (const unsigned int j, double new_pos[]);
    void destroy_segment (const int iconn_index);
    int connection_search (const double pos[]);
    int connection_search_closest (const double pos[]);
    void create_bifurcation (const int iconn_index, Point new_point);

    void update_points (const unsigned int index);
    void update_segments (const int s_index, const unsigned int p_index);
    void destroy_offspring (const int s_index);

    void get_feasible_point (Point *p, const double radius);
    void generate_point_inside_perfusion_area (Point *p, const double radius);
    void generate_point_inside_perfusion_area (double pos[], const double radius);

    bool is_inside_perfusion_area (const Point *p, const double radius);
    bool is_inside_perfusion_area (const double pos[], const double radius);

    void print_points ();
    void print_segments ();
    void write_to_vtk ();
    void draw_perfusion_area (const double radius);

private:
    void calc_middle_segment (double pos[], Segment *s);
    double calc_dthreashold (const double radius, const int num_term);
    bool has_collision (const double pos[], const int iconn);
    void calc_bifurcation_ratio (Segment *s);
    bool collision_detect (Point p1, Point p2, Point p3, Point p4);
};

double calc_size_segment (const Point *p1, const Point *p2);
double calc_poisseulle (const double Q, const double p, const double l);
void generate_point_inside_circle (Point *p, const double radius);
void print_point (const Point p);
double calc_euclidean_dist (const double a[], const double b[]);

#endif