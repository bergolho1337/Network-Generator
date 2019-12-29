#ifndef RANDOM_H
#define RANDOM_H

#include <iostream>
#include <fstream>
#include <string>
#include <cstdio>
#include <cstdlib>
#include <cmath>

#include "../config/config.h"
#include "dSFMT.h"

const uint32_t RANDOM_ARRAY_SIZE = 8000000;

class Random_Generator
{
public:
    uint32_t counter;
    double seed;
    double *array;

public:
	Random_Generator (const uint32_t seed);
    ~Random_Generator ();

    double get_value ();
    void generate ();
	void print ();
};

#endif
