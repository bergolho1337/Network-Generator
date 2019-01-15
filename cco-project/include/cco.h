#ifndef CCO_H
#define CCO_H

#include <cmath>
#include <ctime>
#include <vector>

#include "options.h"

using namespace std;

struct Point
{
    double x, y, z;
}typedef Point;

struct Segment
{
    int type;

    double beta_l;
    double beta_r;
    double Q;

    Point *p1;
    Point *p2;

    Segment *left;
    Segment *right;
    Segment *parent;
}typedef Segment;

class CCO_Network
{
private:
    int num_terminals;

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
};

void generate_point_inside_circle (Point *p, const double radius);

/*
#include "graph.h"

bool is_inside_circle (const double x, const double y, const double r);
void make_root (Graph *the_network, const double Q_perf, const double p_perf);
void generate_terminals (Graph *the_network, const int N_term);
void grow_cco_tree (Graph *the_network, User_Options *options);
void add_terminal (Graph *the_network);
*/

#endif