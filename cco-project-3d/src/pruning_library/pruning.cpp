#include "pruning.h"

SET_PRUNING_FUNCTION (hyperbolic_tangent)
{
    // Get cost function parameters
    bool ret;
    double A;
    ret = get_parameter_value_from_map(config->params,"A",&A);
    double B;
    ret = get_parameter_value_from_map(config->params,"B",&B);
    double C;
    ret = get_parameter_value_from_map(config->params,"C",&C);
    double D;
    ret = get_parameter_value_from_map(config->params,"D",&D);

    //return 50.0*tanh(-0.25*level+3.0) + 50.0;
    return A*tanh(B*level+C) + D;
}

SET_PRUNING_FUNCTION (exponential)
{
    // Get cost function parameters
    bool ret;
    double A;
    ret = get_parameter_value_from_map(config->params,"A",&A);
    double B;
    ret = get_parameter_value_from_map(config->params,"B",&B);
    double C;
    ret = get_parameter_value_from_map(config->params,"C",&C);
    
    //return 100.0*expf(-0.2*level) + 0.0;
    return A*expf(B*level) + C;
}