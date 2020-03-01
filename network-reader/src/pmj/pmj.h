#ifndef PMJ_H_
#define PMJ_H_

#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cassert>
#include <cstring>
#include <string>
#include <vector>
#include <queue>
#include <map>
#include <algorithm>
#include <fstream>

using namespace std;

// =============================================================================================================
// =============================================================================================================
// Structure for a Purkinje-Muscle-Junction
class PMJ
{
public:
    PMJ (const uint32_t index, const double x, const double y, const double z);
    void print ();

public:
	uint32_t index;			// Identifier 
	double pos[3];          // Coordinates
};

// =============================================================================================================
// =============================================================================================================
// Auxiliary functions
void read_pmjs_from_file (vector<PMJ> &pmjs);
// =============================================================================================================

#endif
