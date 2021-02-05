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
#include <stack>
#include <algorithm>

#include "../point/point.h"
#include "../segment/segment.h"
#include "../options/user_options.h"
#include "../utils/utils.h"
#include "../utils/stop_watch.h"
#include "../test/test.h"

#include "../local_optimization_library/local_optimization.h"
#include "../cost_function_library/custom_function.h"

#include "cloud.h"
#include "pmj.h"
#include "helper.h"

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
    Cloud *cloud_data;

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
    uint32_t cur_pmj_package;
    double epsilon_2ms;
    double epsilon_5ms;
    double rmse;
    double rrmse;
    double max_lat_error;
    double min_max_aprox_lat[2];
    double min_max_ref_lat[2];
    double min_term_lat;
    double max_term_lat;

public:
    CCO_Network ();
    CCO_Network (User_Options *options);
    ~CCO_Network ();
    void grow_tree (User_Options *options);
    Segment* build_segment (LocalOptimizationConfig *local_opt_config, const uint32_t index, Point *new_term);
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
    void calc_electric_error ();
    void grow_tree_using_cloud_points (User_Options *options);
    void make_root_using_cloud_points ();
    void make_root_using_initial_network ();
    void rescale_root (Segment *iroot);
    void rescale_until_root (Segment *ipar, Segment *ipar_left, Segment *ipar_right);
    bool generate_terminal_using_cloud_points (CostFunctionConfig *cost_function_config, LocalOptimizationConfig *local_opt_config);
    bool generate_terminal_using_point (CostFunctionConfig *cost_function_config, LocalOptimizationConfig *local_opt_config, Point *pmj_point, const bool evaluate);
    bool generate_terminal (CostFunctionConfig *cost_function_config, LocalOptimizationConfig *local_opt_config);
    bool evaluate_pmj_local_activation_time (Segment *inew, Point *pmj_point, CostFunctionConfig *cost_function_config);
    bool attempt_pmj_connection (CostFunctionConfig *cost_function_config, LocalOptimizationConfig *local_opt_config);
    bool attempt_connect_using_region_radius (Point *pmj_point, CostFunctionConfig *cost_function_config);
    bool attempt_connect_using_inverse (CostFunctionConfig *cost_function_config, LocalOptimizationConfig *local_opt_config, std::vector<Segment*> feasible_segments, Point *pmj_point);
    bool connect_using_local_backtracking (CostFunctionConfig *cost_function_config, LocalOptimizationConfig *local_opt_config, Point *pmj_point);
    bool connect_remaining_active_pmjs (CostFunctionConfig *cost_function_config, LocalOptimizationConfig *local_opt_config);
    bool connect_remaining_inactive_pmjs (CostFunctionConfig *cost_function_config, LocalOptimizationConfig *local_opt_config);
    void prune_segment (Segment *inew);
    bool force_pmj_connection (CostFunctionConfig *cost_function_config, LocalOptimizationConfig *local_opt_config);
    bool force_connection_to_closest_segment (CostFunctionConfig *cost_function_config, LocalOptimizationConfig *local_opt_config, Point *pmj);
    bool attempt_connect_using_local_backtracking (CostFunctionConfig *cost_function_config, LocalOptimizationConfig *local_opt_config, Point *pmj_point);
    bool adjust_terminal_diameter (Segment *inew, Point *pmj_point);
    bool try_connect_pmj (CostFunctionConfig *cost_function_config, LocalOptimizationConfig *local_opt_config, Point *pmj_point);
    bool connect_active_pmjs (CostFunctionConfig *cost_function_config, LocalOptimizationConfig *local_opt_config);
    bool check_active_pmjs_connection ();
    void update_min_max_terminal_lat ();

};

#endif