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

SET_PRUNING_FUNCTION (length)
{
    // Get cost function parameters
    bool ret;
    double length_limit;
    ret = get_parameter_value_from_map(config->params,"length_limit",&length_limit);
    
    if (length < length_limit)
        return 100.0;
    else
        return 0.0;
}

SET_PRUNING_FUNCTION (hyperbolic_tangent_with_length)
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
    double length_limit;
    ret = get_parameter_value_from_map(config->params,"length_limit",&length_limit);

    if (length < length_limit)
        return 100.0;
    else
        return A*tanh(B*level+C) + D;
}