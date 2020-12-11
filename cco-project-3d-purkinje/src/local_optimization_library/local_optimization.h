#ifndef LOCAL_OPTIMIZATION_H
#define LOCAL_OPTIMIZATION_H

#include <cstdio>
#include <cstdlib>
#include <cmath>

#include <string>
#include <vector>

class CCO_Network;
class Segment;
class Point;

// =========================================================================================================
// CONSTANTS AND MACROS
static const uint32_t NE = 5;              // Number of points inside the bifurcation plane (Default = 5)
// =========================================================================================================

// Abstract class - Interface
class LocalOptimization
{
public:
    virtual void optimize (Segment *iconn, Segment *ibiff, Segment *inew,\
                        std::vector<Point*> &test_positions) = 0;
    void move_bifurcation_location (Segment *iconn, Segment *ibiff, Segment *inew, Point *p);
    void move_bifurcation_location (Segment *iconn, Segment *ibiff, Segment *inew, const double pos[]);
    void free_test_positions (std::vector<Point*> &test_positions);
};

// Types of cost function
class DefaultOptimization : public LocalOptimization
{
public:
    void optimize (Segment *iconn, Segment *ibiff, Segment *inew,\
                    std::vector<Point*> &test_positions);
};

class RafaelOptimization : public LocalOptimization
{
public:
    void optimize (Segment *iconn, Segment *ibiff, Segment *inew,\
                    std::vector<Point*> &test_positions);
    bool is_vertex (const uint32_t i, const uint32_t j);
};

#endif