#include "utils.h"

double euclidean_norm (const double x1, const double y1, const double z1,\
                    const double x2, const double y2, const double z2)
{
    return sqrt( pow(x2-x1,2) + pow(y2-y1,2) + pow(z2-z1,2) );
}

double generate_random_number ()
{
    double number = (double)rand() / (double)RAND_MAX;
    double sign = rand() % 2;

    return (sign) ? number : -number;
}

void generate_point_inside_perfusion_area (double pos[], const double radius)
{
    double teta = generate_random_number()*2.0*M_PI;
    double r = fabs(generate_random_number())*radius;

    // Center of perfusion area = (0,-radius,0)
    pos[0] = 0 + r*cos(teta);
    pos[1] = -radius + r*sin(teta);
    pos[2] = 0;
}

double calc_dthreashold (const double radius, const int num_terminals)
{
    return sqrt( (double)(M_PI*radius*radius) / (double)(num_terminals) );
}

double calc_dproj (struct segment_node *s, const double pos[])
{
    struct point *proximal = s->value->src->value;
    struct point *distal = s->value->dest->value;

    double length = euclidean_norm(proximal->x,proximal->y,proximal->z,\
                                   distal->x,distal->y,distal->z);

    double dot_product = (proximal->x - distal->x)*(pos[0] - distal->x) +\
                         (proximal->y - distal->y)*(pos[1] - distal->y) +\
                         (proximal->z - distal->z)*(pos[2] - distal->z);

    return dot_product * pow(length,-2.0);
}

double calc_dortho (struct segment_node *s, const double pos[])
{
    struct point *proximal = s->value->src->value;
    struct point *distal = s->value->dest->value;

    double length = euclidean_norm(proximal->x,proximal->y,proximal->z,\
                                   distal->x,distal->y,distal->z);
    double length_term = euclidean_norm(pos[0],pos[1],pos[2],\
                                distal->x,distal->y,distal->z);

    double dot_product = (proximal->x - distal->x)*(distal->x - pos[0]) +\
                         (proximal->y - distal->y)*(distal->y - pos[1]) +\
                         (proximal->z - distal->z)*(distal->z - pos[2]);

    return sqrt( powf(length_term * length, 2.0) - powf(dot_product,2.0) ) * powf(length,-1.0);
}

double calc_dend (struct segment_node *s, const double pos[])
{
    struct point *proximal = s->value->src->value;
    struct point *distal = s->value->dest->value;

    double d_dist = euclidean_norm(pos[0],pos[1],pos[2],\
                                       distal->x,distal->y,distal->z);

    double d_prox = euclidean_norm(pos[0],pos[1],pos[2],\
                                    proximal->x,proximal->y,proximal->z);

    return std::min(d_dist,d_prox);
}

/*
void draw_perfusion_area (struct cco_network *the_network)
{
    int K_term = the_network->num_terminals;
    int N_term = the_network->N_term;
    double r_perf = the_network->r_perf;
    double A_perf = the_network->A_perf;

    // Increase support domain
    double A_supp = (double)((K_term + 1) * A_perf) / (double)N_term;
    double r_supp = sqrt(A_supp/M_PI);

    // Create a circle
    vtkSmartPointer<vtkRegularPolygonSource> polygonSource =
      vtkSmartPointer<vtkRegularPolygonSource>::New();

    polygonSource->GeneratePolygonOff(); // Outline of the circle
    polygonSource->SetNumberOfSides(50);
    polygonSource->SetRadius(r_supp);
    polygonSource->SetCenter(0,-r_supp,0);

    vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
    writer->SetFileName("output/perfusion_area.vtp");
    writer->SetInputConnection(polygonSource->GetOutputPort());
    writer->Write();
}
*/

// Return true if we detect a collision or if the points are too close
// TODO: Double-check this function again ...
bool collision_detection (const double x1, const double y1, const double z1,\     
                          const double x2, const double y2, const double z2,\
                          const double x3, const double y3, const double z3,\
                          const double x4, const double y4, const double z4)
{
    //printf("p1 = (%g,%g,%g)\n",x1,y1,z1);
    //printf("p2 = (%g,%g,%g)\n",x2,y2,z2);
    //printf("p3 = (%g,%g,%g)\n",x3,y3,z3);
    //printf("p4 = (%g,%g,%g)\n",x4,y4,z4);

    bool ret = true;
    double p13[3], p43[3], p21[3];

    p13[0] = x1 - x3;
    p13[1] = y1 - y3;
    p13[2] = z1 - z3;

    p43[0] = x4 - x3;
    p43[1] = y4 - y3;
    p43[2] = z4 - z3;

    p21[0] = x2 - x1;
    p21[1] = y2 - y1;
    p21[2] = z2 - z1;
    
    if (check_size(p43))
        ret = false;
    if (check_size(p21))
        ret = false;
    
    // Calculate dot products
    double dot_product_1343 = calc_dot_product(p13,p43);
    double dot_product_4321 = calc_dot_product(p43,p21);
    double dot_product_1321 = calc_dot_product(p13,p21);
    double dot_product_4343 = calc_dot_product(p43,p43);
    double dot_product_2121 = calc_dot_product(p21,p21);

    double denominator = dot_product_2121 * dot_product_4343 - dot_product_4321 * dot_product_4321;
    if (fabs(denominator) < EPSILON)
        ret = false;
    
    double numerator = dot_product_1343 * dot_product_4321 - dot_product_1321 * dot_product_4343;

    double mua = numerator / denominator;
    double mub = (dot_product_1343 + dot_product_4321 * mua) / dot_product_4343;

    double pa[3], pb[3];
    pa[0] = x1 + mua * p21[0];
    pa[1] = y1 + mua * p21[1];
    pa[2] = z1 + mua * p21[2];

    pb[0] = x3 + mub * p43[0];
    pb[1] = y3 + mub * p43[1];
    pb[2] = z3 + mub * p43[2];

    //printf("Collision = %d\n\n",ret);
    return ret;
}

void build_unitary_vector (double u[], const double x1, const double y1, const double z1,\
                                       const double x2, const double y2, const double z2)
{
    double norm = euclidean_norm(x1,y1,z1,x2,y2,z2);

    u[0] = (x2 - x1) / norm;
    u[1] = (y2 - y1) / norm;
    u[2] = (z2 - z1) / norm;
}

double calc_angle_between_vectors (const double u[], const double v[])
{
    double dot_product = u[0]*v[0] + u[1]*v[1] + u[2]*v[2];

    double angle_radians = acos(dot_product);

    // Return the angle in degrees
    return angle_radians * 180.0 / M_PI;

}

// Calculate the perfusion radius using the formula of volume of a sphere
double calc_perfusion_radius (const double V)
{
    return powf(((3.0*V)/(4.0 * M_PI)), (1.0/3.0));
}

// Calculate the flux over the terminals. The influx is equally divided for each terminal
double calc_flux_terminals (const double Q_perf, const int N_term)
{
    return Q_perf / N_term;
}

// Calculate the dot product between two vectors
double calc_dot_product (const double u[], const double v[])
{
    return (u[0] * v[0]) + (u[1] * v[1]) + (u[2] * v[2]); 
}

void generate_random_array_using_mersenne_twister (std::vector<double> &the_array)
{
    dsfmt_t dsfmt;
    dsfmt_init_gen_rand(&dsfmt, RAND_SEED);
    
    for (uint32_t i = 0; i < RAND_ARRAY_SIZE; i++)
    {
        double value = dsfmt_genrand_close_open(&dsfmt);
        the_array.push_back(value);
    }
}

bool check_size (const double p[])
{
    if (fabs(p[0]) < EPSILON && fabs(p[1]) < EPSILON && fabs(p[2]) < EPSILON)
        return true;
    else
        return false;
}

void print_terminal_activation_time (struct cco_network *the_network,\
                            const double c, const double cm, const double rc, const double rm)
{
    struct segment_list *s_list = the_network->segment_list;
    struct segment_node *tmp = s_list->list_nodes;

    while (tmp != NULL)
    {
        if (is_terminal(tmp))
        {
            double at = calc_terminal_activation_time(tmp,c,cm,rc,rm);

            //printf("Terminal %d -- AT = %g ms\n",tmp->id,at);
        }
        tmp = tmp->next;
    }
}

void write_to_vtk (struct cco_network *the_network)
{
    uint32_t num_points = the_network->point_list->num_nodes;
    struct point_list *p_list = the_network->point_list;
    struct point_node *p_tmp = p_list->list_nodes;
    uint32_t num_segments = the_network->segment_list->num_nodes;
    struct segment_list *s_list = the_network->segment_list;
    struct segment_node *s_tmp = s_list->list_nodes;

    FILE *file = fopen("output/cco_tree_cm.vtk","w+");
    FILE *file_converted = fopen("output/cco_tree_um.vtk","w+");

    // Write the header
    // File number 1 {cm}
    fprintf(file,"# vtk DataFile Version 3.0\n");
    fprintf(file,"Tree\n");
    fprintf(file,"ASCII\n");
    fprintf(file,"DATASET POLYDATA\n");
    fprintf(file,"POINTS %u float\n",num_points);
    // File number 2 {um}
    fprintf(file_converted,"# vtk DataFile Version 3.0\n");
    fprintf(file_converted,"Tree converted\n");
    fprintf(file_converted,"ASCII\n");
    fprintf(file_converted,"DATASET POLYDATA\n");
    fprintf(file_converted,"POINTS %u float\n",num_points);
    while (p_tmp != NULL)
    {
        fprintf(file,"%g %g %g\n",p_tmp->value->x,p_tmp->value->y,p_tmp->value->z);
        fprintf(file_converted,"%g %g %g\n",p_tmp->value->x*CM_TO_UM,p_tmp->value->y*CM_TO_UM,p_tmp->value->z*CM_TO_UM);
        p_tmp = p_tmp->next;
    }

    fprintf(file,"LINES %u %u\n",num_segments,num_segments*3);
    fprintf(file_converted,"LINES %u %u\n",num_segments,num_segments*3);
    while (s_tmp != NULL)
    {
        fprintf(file,"2 %u %u\n",s_tmp->value->src->id,s_tmp->value->dest->id);
        fprintf(file_converted,"2 %u %u\n",s_tmp->value->src->id,s_tmp->value->dest->id);
        s_tmp = s_tmp->next;
    }
    fprintf(file,"CELL_DATA %u\n",num_segments);
    fprintf(file,"SCALARS radius float\n");
    fprintf(file,"LOOKUP_TABLE default\n");
    fprintf(file_converted,"CELL_DATA %u\n",num_segments);
    fprintf(file_converted,"SCALARS radius float\n");
    fprintf(file_converted,"LOOKUP_TABLE default\n");
    s_tmp = s_list->list_nodes;
    while (s_tmp != NULL)
    {
        fprintf(file,"%g\n",s_tmp->value->radius * 1000.0);
        fprintf(file_converted,"%g\n",s_tmp->value->radius * 1000.0);
        s_tmp = s_tmp->next;
    }
    fclose(file);
}
