#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>

#include <vector>
#include <string>
#include <list>
#include <algorithm>

// ============================================================================================================================================
// CONSTANTS AND MACROS 
// ============================================================================================================================================
#define UM_TO_CM 0.0001
#define CM_TO_UM 1.0E+05
#define MM_TO_UM 1.0E+03
#define M_TO_UM 1.0E+06
#define M_TO_MM 1.0E+03
#define MS_TO_S 0.001
#define S_TO_MS 1000.0
#define CM_S_TO_M_S 0.01
#define M_S_TO_UM_MS 1000.0
#define MS_TO_US 1000.0
#define CM_TO_MM 10.0
#define CM_TO_M 0.01

static const double ETA = 3.6e-03;                  // Blood viscosity 
static const uint32_t NTOSS = 10;                   // Number of tosses for a new terminal
static const uint32_t NCONN = 20;                   // Number of segments to test for connection
static const uint32_t PRUNING_PASSES = 1;           // Number of times the pruning procedure will be called
static const uint32_t MAX_PMJ_PACKAGE_TRY = 1000;   // Maximum number of times we will try to connect the points inside the PMJ package
static const double FACTOR = 0.95;                  // Reduction factor for the distance criterion
static const double D_THREASH_LIMIT = 1.0E-10;      // Limit for the d_threash
static const double PMJ_LOOSE_RATE = 1.2;           // Loose rate for the PMJ {increase by 20%} (forcing connection)
static const double MIN_CV_THREASHOLD = 1000.0;     // Minimum threashold for the conduction velocity {um/ms} (diameter calibration)
static const double MAX_CV_THREASHOLD = 4000.0;     // Maximum threashold for the conduction velocity {um/ms} (diameter calibration)
// ============================================================================================================================================

#endif
