//
// Created by bergolho on 12/02/19.
//

#ifndef CCO_H
#define CCO_H

#include <cstdio>
#include <cstdlib>
#include <cstdbool>
#include <cstdint>
#include <cmath>

#include <vector>
#include <string>
#include <algorithm>

#include "../point/point.h"
#include "../segment/segment.h"
#include "../options/user_options.h"
#include "../utils/utils.h"
#include "../utils/stop_watch.h"
#include "../test/test.h"

#include "../local_optimization_library/local_optimization.h"
#include "../cost_function_library/cost_function.h"

#include "cloud_points.h"
#include "pmj.h"

// CONSTANTS AND MACROS 
// ============================================================================================================================================
static const double ETA = 3.6e-03;                  // Blood viscosity 
static const uint32_t NTOSS = 10;                   // Number of tosses for a new terminal
static const uint32_t NCONN = 20;                   // Number of segments to test for connection
static const double FACTOR = 0.95;                  // Reduction factor for the distance criterion
static const double D_THREASH_LIMIT = 1.0e-05;      // Limit for the d_threash
static const uint32_t PRUNING_PASSES = 1;           // Number of times the pruning procedure will be called
static const uint32_t PMJ_LOOSE_THREASHOLD = 5.0;   // Threashold for loosening the LAT error tolerance and forcing the PMJ connection {ms}
static const double PMJ_LOOSE_RATE = 1.2;           // Loose rate for the PMJ {increase by 20%} (forcing connection)
static const double CV_THREASHOLD = 1000.0;         // Threashold for the conduction velocity {um/ms} (diameter calibration)
// ============================================================================================================================================

class CCO_Network
{
public:
    // Constant declarations
    static const double Q_PERF;                  
    static const double P_PERF;                     
    static const double P_TERM;                      
    static const double V_PERF;
    double R_PERF;                    
    double Q_TERM;
    uint32_t CUR_MAX_PMJ_INDEX;

    // Variable declarations
    uint32_t num_terminals;
    uint32_t seed;
    uint32_t max_rand_offset;
    uint32_t N_term;

    bool using_only_murray_law;
    double start_radius;
    double gamma;
    double lat_offset;

    double root_pos[3];

    std::vector<Point*> point_list;
    std::vector<Segment*> segment_list;

    bool using_cloud_points;
    std::string cloud_points_filename;
    std::vector<bool> cloud_points_connected;
    std::vector<Point*> cloud_points;

    bool using_obstacle;
    //std::vector<Face> obstacles;
    
    bool using_pmj_location;
    PMJ *pmj_data;

    std::string cost_function_name;
    CostFunction *cost_fn;

    bool using_local_optimization;
    std::string local_optimization_function_name;
    LocalOptimization *local_opt_fn;

    //bool using_pruning;
    //Pruning *prun;

    std::string output_dir;
    FILE *log_file;
private:
    uint32_t cur_rand_index;
    double epsilon_2ms;
    double epsilon_5ms;
    double rmse;
    double rrmse;
    double max_lat_error;
    double min_max_aprox_lat[2];
    double min_max_ref_lat[2];

public:
    CCO_Network ();
    CCO_Network (User_Options *options);
    ~CCO_Network ();
    void grow_tree (User_Options *options);
    Segment* build_segment (LocalOptimizationConfig *local_opt_config, const uint32_t index, Point *new_term);
    Segment* get_terminal (const double pos[]);
    Segment* get_terminal_2 (const double pos[]);
    double calc_terminal_local_activation_time (Segment *term);
    CCO_Network* copy ();
    CCO_Network* concatenate (CCO_Network *input);
    void rescale_tree (Segment *ibiff, Segment *iconn, Segment *inew);
    void recalculate_radius ();
    void recalculate_length ();
    void restore_state_tree (Segment *iconn);
    void link_segments (Segment *term_1, Segment *term_2);
    void adjust_radius ();
    void adjust_radius_2 ();
    void print ();
    void print_point_list ();
    void print_segment_list ();
    void print_network_info ();
    void write_to_vtk_iteration ();
private:
    void set_parameters (User_Options *options);
    void set_cost_function ();
    void set_save_network ();
    void set_local_optimization_function_name ();
    void get_segment_length (std::vector<double> &segments);
    void get_bifurcation_angles(std::vector<double> &angles);
    void calc_electric_error ();
    void read_cloud_points ();
    void read_pmj_locations (PMJConfig *config);
    void grow_tree_using_cloud_points (User_Options *options);
    void make_root_using_cloud_points ();
    void make_root_using_initial_network ();
    void rescale_root (Segment *iroot);
    void rescale_until_root (Segment *ipar, Segment *ipar_left, Segment *ipar_right);
    void update_segment_pointers (Segment *iconn, Segment *ibiff, Segment *inew, const bool is_restore);
    void update_ndist (Segment *ibiff, const bool is_restore); 
    bool generate_terminal_using_cloud_points (CostFunctionConfig *cost_function_config, LocalOptimizationConfig *local_opt_config);
    bool generate_terminal_using_pmj_locations (CostFunctionConfig *cost_function_config, LocalOptimizationConfig *local_opt_config, Point *pmj_point, const bool evaluate);
    bool connection_search (Point *p, const double d_threash);
    bool distance_criterion (Segment *s, const double pos[], const double d_threash);
    bool check_null_segments ();
    bool check_bifurcation_rule ();
    bool check_collisions_and_fill_feasible_segments (Point *p, std::vector<Segment*> &feasible_segments);
    bool has_collision (Segment *s, Point *p);
    bool has_collision (Segment *iconn, Segment *ibiff, Segment *inew);
    bool evaluate_pmj_local_activation_time (Segment *inew, Point *pmj_point, CostFunctionConfig *cost_function_config);
    bool attempt_pmj_connection (CostFunctionConfig *cost_function_config, LocalOptimizationConfig *local_opt_config);
    bool attempt_connect_using_region_radius (Point *pmj_point, CostFunctionConfig *cost_function_config);
    bool attempt_connect_using_max_lat (CostFunctionConfig *cost_function_config, LocalOptimizationConfig *local_opt_config, Point *pmj_point);
    bool force_pmj_connection (CostFunctionConfig *cost_function_config, LocalOptimizationConfig *local_opt_config, Point *pmj_point);
    bool force_connection_to_closest_segment (CostFunctionConfig *cost_function_config, LocalOptimizationConfig *local_opt_config, Point *pmj_point);
    bool connect_remaining_active_pmjs (CostFunctionConfig *cost_function_config, LocalOptimizationConfig *local_opt_config);
    bool connect_remaining_inactive_pmjs (CostFunctionConfig *cost_function_config, LocalOptimizationConfig *local_opt_config);
    Point* generate_bifurcation_node (Segment *iconn, LocalOptimizationConfig *local_opt_config);
    Point* generate_terminal_node (Point *p);
    Point* search_point (const uint32_t index);
    void eliminate_point_from_list (Point *p);
    void eliminate_segment_from_list (Segment *s);
    void prune_segment (Segment *inew);
    void order_point_list ();
    void order_segment_list ();
    uint32_t sort_point_from_cloud (Point *p);
};

#endif