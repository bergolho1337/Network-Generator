// Author: Lucas Berg

#include <iostream>
#include <vector>
#include <random>
#include <algorithm>

#include <cstdio>
#include <cstdlib>
#include <cmath>

using namespace std;

// ====================================================================================================================
// CONSTANTS AND MACROS
#define PRINT_LINE "--------------------------------------------------------------------------------------------"

const double MEAN = 1.0;
const double STD_DEV = 0.1;
const uint32_t NPOINTS = 1000; 
// ====================================================================================================================

void generate_random_numbers (vector<double> &arr)
{
    default_random_engine generator;
    normal_distribution<double> distribution(MEAN,STD_DEV);

    for (uint32_t i = 0; i < NPOINTS; i++)
    {
        double number = distribution(generator);
        
        arr.push_back(number);
    }
}

void print_array (const vector<double> arr)
{
    for (uint32_t i = 0; i < arr.size(); i++)
        printf("Array %u = %g -- %d\n",i,arr[i],(int)arr[i]);
}

void print_histogram (const vector<double> arr)
{
    // Create an array with 0 entries
    uint32_t p[10] = {};

    for (uint32_t i = 0; i < arr.size(); i++)
    {
        if (arr[i] >= 0.0 && arr[i] <= 10.0)
            p[(int)arr[i]]++;
    }

    for (uint32_t i = 0; i < 10; i++)
    {
        printf("%02u-%02u: %s\n",i,i+1,string(p[i]*100/NPOINTS,'*').c_str());
    }
        
}

double get_number_from_normal_distribution (vector<double> arr)
{
    uint32_t index = rand() % arr.size();

    return arr[index];
}

int main (int argc, char *argv[])
{
    if (argc-1 != 0)
    {
        printf("%s\n",PRINT_LINE);
        printf("Usage:> %s\n",argv[0]);
        printf("%s\n",PRINT_LINE);
        exit(EXIT_FAILURE);   
    }

    vector<double> arr;

    generate_random_numbers(arr);
    
    print_array(arr);
    
    print_histogram(arr);

    srand(time(NULL));
    double start_radius = 1.0;
    double k = 3.0;
    double N = get_number_from_normal_distribution(arr);
    double r_left = N * pow( 0.5*pow(start_radius,k) , 1/k );
    double r_right = pow ( pow(start_radius,k) - pow(r_left,k) , 1/k );

    printf("r_0 = %g\n",start_radius);
    printf("r_left = %g\n",r_left);
    printf("r_right = %g\n",r_right);

    return 0;
}