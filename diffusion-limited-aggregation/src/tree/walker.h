//
// Created by bergolho on 26/01/20.
//

#ifndef WALKER_H
#define WALKER_H

#include <cstdio>
#include <cstdlib>
#include <cstdbool>
#include <cstdint>
#include <cmath>
#include <vector>
#include <string>

const double WIDTH = 600.0;
const double HEIGHT = 600.0;
const double RADIUS = 6.0;

class Walker
{
public:
    double pos[3];
    bool stuck;
    double radius;
public:
    Walker ();
    Walker (const double x, const double y, const double z);
    void walk ();
    uint32_t is_stuck (std::vector<Walker*> the_tree);
    void print ();
};

void sort_random_points (double pos[]);
bool is_inside (const double x, const double y);
double calculate_distance (const double p1[], const double p2[]);
double generate_random_number ();
void print_progress_bar (const uint32_t cur_iter, const uint32_t max_iter);

#endif