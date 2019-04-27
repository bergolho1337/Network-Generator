//
// Created by bergolho on 28/02/19.
//

#ifndef TEST_H
#define TEST_H

#include <cstdio>
#include <cstdlib>
#include <cstdbool>
#include <cstdint>
#include <cmath>

#include "../cco/cco.h"

static const double TOLERANCE = 1.0e-05;

// Test functions
void test1 (struct cco_network *the_network);
void test2 (struct cco_network *the_network); 
void test3 (struct cco_network *the_network);

// Unitary test
void check_bifurcation_rule (struct cco_network *the_network);

#endif