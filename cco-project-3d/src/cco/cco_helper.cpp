#include "cco_helper.h"

void calc_middle_point_segment (struct segment_node *s, double pos[])
{
    struct point *src = s->value->src->value;
    struct point *dest = s->value->dest->value;

    pos[0] = (src->x + dest->x) / 2.0;
    pos[1] = (src->y + dest->y) / 2.0;
    pos[2] = (src->z + dest->z) / 2.0;

    printf("prox = (%g,%g,%g) -- distal = (%g,%g,%g) -- meio = (%g,%g,%g)\n",src->x,src->y,src->z,dest->x,dest->y,dest->z,pos[0],pos[1],pos[2]);
}

void calc_relative_resistance_subtree (struct segment_node *ibiff, struct segment_node *iconn, struct segment_node *inew)
{
    calc_relative_resistance_term(ibiff);
    double R = ibiff->value->resistance;
    double R_left = pow(iconn->value->beta,4) / iconn->value->resistance;
    double R_right = pow(inew->value->beta,4) / inew->value->resistance;

    ibiff->value->resistance = R + pow( R_left + R_right , -1.0 );
}

void calc_relative_resistance_term (struct segment_node *iterm)
{
    struct point *src = iterm->value->src->value;
    struct point *dest = iterm->value->dest->value;
    double length = euclidean_norm(src->x,src->y,src->z,dest->x,dest->y,dest->z);

    iterm->value->resistance = 8.0 * ETA * length / M_PI;
}

void calc_pressure_drop_term (struct segment_node *iterm, const double Q_term)
{
    iterm->value->delta_p = iterm->value->resistance * Q_term;
}

void calc_pressure_drop_subtree (struct segment_node *iconn, const double Q_term)
{
    double Q_iconn = iconn->value->ndist * Q_term;

    iconn->value->delta_p = iconn->value->resistance * Q_iconn;
}

void calc_radius_term (struct segment_node *iterm, const double Q_term, const double delta_p)
{
    double R = iterm->value->resistance;
    double Q = iterm->value->ndist * Q_term;

    iterm->value->radius = pow( R * Q / delta_p , 0.25 );
}

double calc_radius_ratio (struct segment_node *iconn, struct segment_node *inew, const double Q_term)
{
    double Q_inew = inew->value->ndist*Q_term;
    double R_inew = inew->value->resistance;

    double Q_iconn = iconn->value->ndist*Q_term;
    double R_iconn = iconn->value->resistance;

    double ratio = pow( (Q_iconn * R_iconn) / (Q_inew * R_inew), 0.25 );

    return ratio;
}

double calc_bifurcation_ratio (const double radius_ratio, const bool sign)
{
    static const double expoent = -1.0/GAMMA;
    double base = radius_ratio;

    if (sign)
        return pow( 1.0 + ( pow( base, -GAMMA) ) , expoent);
    else
        return pow( 1.0 + ( pow( base, GAMMA) ) , expoent);
}

double calc_segment_volume (struct segment_node *s)
{
    struct point *src = s->value->src->value;
    struct point *dest = s->value->dest->value;

    double l = euclidean_norm(src->x,src->y,src->z,\
                              dest->x,dest->y,dest->z);
    double r = s->value->radius;

    return M_PI * r * r * l;
}

double calc_tree_volume (struct cco_network *the_network)
{
    struct segment_list *s_list = the_network->segment_list;
    struct segment_node *tmp = s_list->list_nodes;

    double total_volume = 0.0;

    while (tmp != NULL)
    {
        total_volume += calc_segment_volume(tmp);

        tmp = tmp->next;
    }

    return total_volume;
}

double calc_assymetric_ratio (struct segment_node *right, struct segment_node *left)
{
    // Calculate assymetric ratio using equation (2.12) from Rafael's thesis
    double r_right = right->value->radius;
    double r_left = left->value->radius;

    double epsilon = std::min(r_right,r_left) / std::max(r_right,r_left);

    return epsilon;
}

// Output will be given in (ms)
double calc_tau_m (const double cm, const double rm)
{
    // Equation (4.16) from Keener's book
    return rm * cm;
}

// Input should be given in: ({cm},{ohm.cm},{kohm.cm^2}
// Output will be given in (cm)
double calc_lambda_m (const double r, const double rc, const double rm)
{
    double d = 2.0 * r;

    // Equation (4.17) from Keener's book
    double num = 1000.0 * rm * d;                    // Here we need to convert the Rm to {ohm.cm^2}
    double den = 4.0 * rc;

    return sqrt(num / den);
}

// Calculate the propagation velocity using the cable equation (6.21) from Kenner's book
double calc_propagation_velocity (const double r,\
                        const double c, const double cm, const double rc, const double rm)
{
    return c * calc_lambda_m(r,rc,rm) / calc_tau_m(cm,rm);
}

double calc_terminal_activation_time (struct segment_node *s,\
                        const double c, const double cm, const double rc, const double rm)
{
    struct segment_node *tmp = s;
    double at = 0.0;

    while (tmp != NULL)
    {
        at += calc_segment_activation_time(tmp,c,cm,rc,rm);

        tmp = tmp->value->parent;
    }

    return at;
}

double calc_segment_activation_time (struct segment_node *s,\
                        const double c, const double cm, const double rc, const double rm)
{
    struct point *src = s->value->src->value;
    struct point *dest = s->value->dest->value;
    double length = euclidean_norm(src->x,src->y,src->z,\
                                  dest->x,dest->y,dest->z);

    double delta_s = length;
    double r = s->value->radius;

    double velocity = calc_propagation_velocity(r,c,cm,rc,rm);

    // Parse every distance to {mm}
    velocity *= CM_TO_MM;
    delta_s *= CM_TO_MM;
    r *= CM_TO_MM;

    //printf("\tPropagation velocity = %g mm/ms \n",velocity);
    //printf("\tDiameter = %g um\n",r*2.0*MM_TO_UM);
    //printf("\tDistance = %g mm \n",delta_s);
    //printf("\tRadius = %g mm\n",r);
    //printf("\tActivation time = %g ms\n\n",delta_s / velocity);

    // The output will be on {mm/ms}
    return delta_s / velocity;
}

// The activation time is already in microseconds
double calc_segment_activation_time_using_level (const double at, struct segment_node *iconn)
{
    double level = calc_segment_level(iconn);

    return pow( at, 1.0/level );
}

double calc_segment_level (struct segment_node *iconn)
{
    struct segment_node *tmp = iconn;

    double level = 1.0;
    while (tmp != NULL)
    {
        level++;

        tmp = tmp->value->parent;
    }

    return level;
}

double calc_segment_custom_function_with_level_penalty (const double eval, struct segment_node *iconn)
{
    double level = calc_segment_level(iconn);

    //return eval / level;
    //return eval / (0.5 * level);
    return pow( eval, 1.0/level );
}

double calc_segment_custom_function (struct segment_node *s, const double beta, const double alpha)
{
    struct point *src = s->value->src->value;
    struct point *dest = s->value->dest->value;

    double l = euclidean_norm(src->x,src->y,src->z,\
                              dest->x,dest->y,dest->z);
    double r = s->value->radius;

    return pow(l,beta) * pow(r,alpha);
}

double calc_custom_function (struct cco_network *the_network, const double beta, const double alpha)
{
    struct segment_list *s_list = the_network->segment_list;
    struct segment_node *tmp = s_list->list_nodes;

    double total_eval = 0.0;

    while (tmp != NULL)
    {
        total_eval += calc_segment_custom_function(tmp,beta,alpha);

        tmp = tmp->next;
    }

    return total_eval;
}

bool has_deviation (struct segment_list *s_list, struct segment_node *inew,\
                    const double new_at, const double limit,\
                    const double c, const double cm, const double rc, const double rm)
{
    struct segment_node *tmp = s_list->list_nodes;

    while (tmp != NULL)
    {
        if (tmp != inew && is_terminal(tmp))
        {
            double at = calc_terminal_activation_time(tmp,c,cm,rc,rm);
            if (fabs(new_at - at) > limit)
                return true;
        }
        tmp = tmp->next;
    }

    return false;
}

bool is_terminal (struct segment_node *s)
{
    if (s->value->left == NULL && s->value->right == NULL)
        return true;
    else
        return false;
}