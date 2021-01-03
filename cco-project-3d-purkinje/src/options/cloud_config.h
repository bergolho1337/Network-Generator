#ifndef CLOUD_CONFIG_H
#define CLOUD_CONFIG_H

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <string>

class Segment;
class Point;

class CloudConfig
{
public:
    std::string filename;
    double rand_offset;
    // Additional parameters ...
public:
    CloudConfig ();
    ~CloudConfig ();
    void print ();
};

#endif