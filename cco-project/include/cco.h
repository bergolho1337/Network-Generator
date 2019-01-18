#ifndef CCO_H
#define CCO_H

#include <cmath>
#include <ctime>
#include <vector>
#include <algorithm>

#include "options.h"

// Constants and Macros
const int ROOT_STICKER = 0;
const int TERMINAL_STICKER = 1;
const int BIFURCATION_STICKER = 2;
const int UNDEFINED_VALUE = -1;
const double ETA = 3.6e-03;         // Viscosity of blood (Pa.s)
const int N_toss = 200;             // Number of tosses for a new terminal

#define PRINT_LINE "------------------------------------------------------------------------------"

using namespace std;

class Point
{
public:
    double x, y, z;

public:
    Point ();
    Point (const double x, const double y, const double z);
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
    double length;

    int index_source;
    int index_destination;
    Point *p1;
    Point *p2;

    Segment *left;
    Segment *right;
    Segment *parent;
public:
    Segment ();
    Segment (Point *p1, Point *p2,\
            int index_source, int index_destination,\
            Segment *left, Segment *right, Segment *parent,\
            const double Q, const double p);
    void set_parent (Segment *new_parent) { this->parent = new_parent; }
    void set_left_offspring (Segment *new_offspring) { this->left = new_offspring; }
    void set_right_offspring (Segment *new_offspring) { this->right = new_offspring; }

    void add_offspring (Segment *new_segment);
    double calc_dproj (const Point p);
    double calc_dortho (const Point p);
    double calc_dend (const Point p);
private:
    int get_segment_type ();
    void calc_bifurcation_ratio ();
};

class CCO_Network
{
private:
    int num_terminals;

    int N_term;
    double Q_perf;
    double p_perf;
    double r_perf;
    
    vector<Point> points;
    vector<Segment> segments;

public:
    CCO_Network ();
    CCO_Network (User_Options *options);
    void grow_tree ();
    void make_root ();
    void generate_new_terminal ();
    void create_bifurcation (const int iconn_index, Point new_point);

    void print_points ();
    void print_segments ();
    void write_to_vtk ();

private:
    void calc_middle_segment (Point *p, const Segment segment);
    double calc_dthreashold (const double radius, const int num_term);
    bool has_collision (const Point p, const unsigned int iconn_index);
};

double calc_size_segment (const Point *p1, const Point *p2);
double calc_poisseulle (const double Q, const double p, const double l);
void generate_point_inside_circle (Point *p, const double radius);
void print_point (const Point p);

/*
#include "graph.h"

bool is_inside_circle (const double x, const double y, const double r);
void make_root (Graph *the_network, const double Q_perf, const double p_perf);
void generate_terminals (Graph *the_network, const int N_term);
void grow_cco_tree (Graph *the_network, User_Options *options);
void add_terminal (Graph *the_network);
*/

#endif