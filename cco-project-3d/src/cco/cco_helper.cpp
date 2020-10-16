#include "cco_helper.h"

void calc_middle_point_segment (struct segment_node *s, double pos[])
{
    struct point *src = s->value->src->value;
    struct point *dest = s->value->dest->value;

    pos[0] = (src->x + dest->x) / 2.0;
    pos[1] = (src->y + dest->y) / 2.0;
    pos[2] = (src->z + dest->z) / 2.0;

    //printf("prox = (%g,%g,%g) -- distal = (%g,%g,%g) -- meio = (%g,%g,%g)\n",src->x,src->y,src->z,dest->x,dest->y,dest->z,pos[0],pos[1],pos[2]);
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

double calc_bifurcation_ratio (const double gamma, const double radius_ratio, const bool sign)
{
    static const double expoent = -1.0/gamma;
    double base = radius_ratio;

    if (sign)
        return pow( 1.0 + ( pow( base, -gamma) ) , expoent);
    else
        return pow( 1.0 + ( pow( base, gamma) ) , expoent);
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

// Cable equation: {m/s}
double calc_propagation_velocity (const double d,\
                        const double G, const double Cf, const double tau_f)
{
    return powf( (G*d)/(4.0*Cf*tau_f) , 0.5 ) * 0.1;
}

double calc_terminal_activation_time (struct segment_node *s,\
                        const double G, const double Cf, const double tau_f)
{
    struct segment_node *tmp = s;
    double at = 0.0;

    at = calc_segment_activation_time(tmp,G,Cf,tau_f);
    while (tmp != NULL)
    {
        at += calc_segment_activation_time(tmp,G,Cf,tau_f);

        tmp = tmp->value->parent;
    }

    return at;
}

double calc_total_terminal_activation_time (struct cco_network *the_network,\
                        const double G, const double Cf, const double tau_f)
{
    struct segment_list *s_list = the_network->segment_list;
    struct segment_node *tmp = s_list->list_nodes;

    double total_activation_time = 0.0;

    while (tmp != NULL)
    {
        if (is_terminal(tmp))
        {
            total_activation_time += calc_terminal_activation_time(tmp,G,Cf,tau_f);
        }

        tmp = tmp->next;
    }

    return total_activation_time;
}

double calc_total_activation_time (struct cco_network *the_network,\
                            const double G, const double Cf, const double tau_f)
{
    struct segment_list *s_list = the_network->segment_list;
    struct segment_node *tmp = s_list->list_nodes;

    double total_activation_time = 0.0;

    while (tmp != NULL)
    {
        total_activation_time += calc_segment_activation_time(tmp,G,Cf,tau_f);

        tmp = tmp->next;
    }

    return total_activation_time;
}

// Activation time: {ms}
double calc_segment_activation_time (struct segment_node *s,\
                        const double G, const double Cf, const double tau_f)
{
    struct point *src = s->value->src->value;
    struct point *dest = s->value->dest->value;
    double length = euclidean_norm(src->x,src->y,src->z,\
                                  dest->x,dest->y,dest->z);

    double delta_s = length;                    // m
    double radius = s->value->radius;           // m 
    double diameter = radius * 2.0;             // m

    diameter *= 1000.0;                       // um
    //diameter *= 1000;                         // um
    //delta_s *= 1;

    double velocity = calc_propagation_velocity(diameter,G,Cf,tau_f);
    
    // The output will be on {s}
    double at = delta_s / velocity * S_TO_MS;
    
    //printf("\tPropagation velocity = %g m/s \n",velocity);
    //printf("\tDiameter = %g um\n",diameter);
    //printf("\tDistance = %g m \n",delta_s);
    //printf("\tActivation time = %g ms\n\n",at);
    //exit(1);

    return at;
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

    double level = 0.0;
    while (tmp != NULL)
    {
        level++;

        tmp = tmp->value->parent;
    }

    return level;
}

double calc_level_threashold (const uint32_t num_terminals)
{
    if (num_terminals < 80)
        return 0.0;
    else
        return 20.0;
}

double calc_segment_custom_function_with_level_penalty (const double eval, struct segment_node *iconn)
{
    double level = calc_segment_level(iconn);

    //return eval;
    return eval / level;
    //return eval / (0.5 * level);
    //return pow( eval, 1.0/level );
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

void calc_unitary_vector (struct segment_node *s, double u[])
{
    struct point *src = s->value->src->value;
    struct point *dest = s->value->dest->value;
    double norm = euclidean_norm(src->x,src->y,src->z,dest->x,dest->y,dest->z);
    if (norm < 1.0e-08)
        norm = 1.0e-08;

    u[0] = (dest->x - src->x) / norm;
    u[1] = (dest->y - src->y) / norm;
    u[2] = (dest->z - src->z) / norm;
}

double calc_segment_size (struct segment_node *s)
{
    struct point *src = s->value->src->value;
    struct point *dest = s->value->dest->value;
    double norm = euclidean_norm(src->x,src->y,src->z,dest->x,dest->y,dest->z);

    return norm;
}

uint32_t calc_closest_pmj_site (struct segment_node *s, std::vector<struct pmj_lat> lat_points)
{
    uint32_t closest_id = lat_points.size()+1;
    double closest_dist = __DBL_MAX__;

    double x1 = s->value->dest->value->x;
    double y1 = s->value->dest->value->y;
    double z1 = s->value->dest->value->z;
    for (uint32_t i = 0; i < lat_points.size(); i++)
    {
        double x2 = lat_points[i].x;
        double y2 = lat_points[i].y;
        double z2 = lat_points[i].z;

        double dist = euclidean_norm(x1,y1,z1,x2,y2,z2);
        if (dist < closest_dist)
        {
            closest_dist = dist;
            closest_id = i;
        }
    }

    return closest_id;
}



bool has_deviation (struct segment_list *s_list, struct segment_node *inew,\
                    const double new_at, const double limit,\
                    const double G, const double Cf, const double tau_f)
{
    struct segment_node *tmp = s_list->list_nodes;

    while (tmp != NULL)
    {
        if (tmp != inew && is_terminal(tmp))
        {
            double at = calc_terminal_activation_time(tmp,G,Cf,tau_f);
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

bool check_angle_restriction (const double angle, const double min_angle, const double max_angle)
{
    return (angle > min_angle && angle < max_angle);
}

void read_initial_network (const char filename[], std::vector<struct cco_point> &the_points, std::vector<struct cco_segment> &the_segments)
{
    FILE *file = fopen(filename,"r");
    if (!file)
    {
        fprintf(stderr,"[cco] ERROR! Reading initial network '%s'\n",filename);
        exit(EXIT_FAILURE);
    }

    uint32_t num_points;
    fscanf(file,"%u",&num_points);
    
    std::vector<struct cco_point> points;
    for (uint32_t i = 0; i < num_points; i++)
    {
        struct cco_point p;
        fscanf(file,"%lf %lf %lf",&p.x,&p.y,&p.z);

        the_points.push_back(p);
    }
    uint32_t num_segments;
    fscanf(file,"%u",&num_segments);

    std::vector<struct cco_segment> segments;
    for (uint32_t i = 0; i < num_segments; i++)
    {
        struct cco_segment s;
        fscanf(file,"%d %d %d",&s.parent_id,&s.left_offspring_id,&s.right_offspring_id);
        fscanf(file,"%d %d",&s.src_id,&s.dest_id);

        the_segments.push_back(s);
    }

    fclose(file);
}