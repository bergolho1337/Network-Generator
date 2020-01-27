//
// Created by bergolho on 26/01/20.
//

#ifndef SEGMENT_H
#define SEGMENT_H

#include <cstdio>
#include <cstdlib>
#include <cstdbool>
#include <cstdint>
#include <vector>

class Segment
{
public:
    uint32_t src;
    uint32_t dest;
public:
    Segment (const uint32_t src, const uint32_t dest);
};


#endif