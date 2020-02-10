#include "utils.h"

struct random_generator* new_random_generator ()
{
    struct random_generator *result = (struct random_generator*)malloc(sizeof(struct random_generator));

    result->counter = 0;
    result->array = NULL;

    return result;
}

void free_random_generator (struct random_generator *the_generator)
{
    if (the_generator->array)
        free(the_generator->array);
    free(the_generator);
}

void generate_random_array (struct random_generator *the_generator)
{
    the_generator->array = (double*)malloc(sizeof(double)*RAND_ARRAY_SIZE);
    
    dsfmt_t dsfmt;
    dsfmt_init_gen_rand(&dsfmt, RAND_SEED);
    for (uint32_t i = 0; i < RAND_ARRAY_SIZE; ++i) 
    {
        the_generator->array[i] = dsfmt_genrand_close_open(&dsfmt);
    }
}

double get_value (struct random_generator *the_generator)
{
    double value = the_generator->array[the_generator->counter];
    
    the_generator->counter++;
    if (the_generator->counter > RAND_ARRAY_SIZE)
        the_generator->counter = 0;
    
    return value;
}

void generate_cloud_points (struct random_generator *the_generator, std::vector<struct point> &cloud_points, const double radius)
{
    printf("\n[utils] Generating default sphere cloud of points (Total number of points = %u)\n",TOTAL_CLOUD_POINTS_SIZE);

    double pos[3];
    for (uint32_t i = 0; i < TOTAL_CLOUD_POINTS_SIZE; i++)
    {
        bool point_is_inside_sphere = false;
        while (!point_is_inside_sphere)
        {
            pos[0] = 2.0 * get_value(the_generator) - 1.0;
            pos[1] = 2.0 * get_value(the_generator) - 1.0;
            pos[2] = 2.0 * get_value(the_generator) - 1.0;

            // Convert to the real domain
            pos[0] *= radius;
            pos[1] *= radius;
            pos[2] *= radius;
            
            if (sqrt(pos[0] * pos[0] + pos[1] * pos[1] + pos[2] * pos[2]) > radius)
                point_is_inside_sphere = false;
            else
                point_is_inside_sphere = true;

        }

        struct point point;
        point.x = pos[0];
        point.y = pos[1];
        point.z = pos[2];

        cloud_points.push_back(point);
        
    }
}

double euclidean_norm (const double x1, const double y1, const double z1,\
                    const double x2, const double y2, const double z2)
{
    return sqrt( pow(x2-x1,2) + pow(y2-y1,2) + pow(z2-z1,2) );
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

void draw_perfusion_volume (const double radius)
{
    // Create a sphere
    vtkSmartPointer<vtkSphereSource> sphere_source = vtkSmartPointer<vtkSphereSource>::New();
    sphere_source->SetCenter(0.0, 0.0, 0.0);
    sphere_source->SetRadius(radius);

    // Make the surface smooth.
    sphere_source->SetPhiResolution(100);
    sphere_source->SetThetaResolution(100);

    vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
    writer->SetFileName("output/perfusion_volume.vtp");
    writer->SetInputConnection(sphere_source->GetOutputPort());
    writer->Write();
}

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

// Calculate the difference between two vectors
void calc_subtract_vector(const double a[], const double b[], double c[])
{
    for (uint32_t i = 0; i < 3; i++)
        c[i] = a[i] - b[i];
}

void calc_normal_vector (const double a[], const double b[], const double c[], double n[])
{
    n[0] = ( c[2] - a[2] ) * ( b[1] - a[1] ) -\
           ( b[2] - a[2] ) * ( c[1] - a[1]);
    n[1] = ( b[2] - a[2] ) * ( c[0] - a[0] ) -\
           ( b[0] - a[0] ) * ( c[2] - a[2]);
    n[2] = ( b[0] - a[0] ) * ( c[1] - a[1] ) -\
           ( b[1] - a[1] ) * ( c[0] - a[0]);
}

void calc_plane_coefficients (const double v1[], const double v2[], const double v3[], double N[], double &D)
{
    // Plane equation: ax + by + cz = d --> N = (a,b,c)
    calc_normal_vector(v1,v2,v3,N);
    D = calc_dot_product(N,v1);
}

bool check_size (const double p[])
{
    if (fabs(p[0]) < EPSILON && fabs(p[1]) < EPSILON && fabs(p[2]) < EPSILON)
        return true;
    else
        return false;
}

bool check_segment_plane_intersection (const double x_prox[], const double x_new[], struct face the_face)
{
    // Calculate the plane equation for the given triangle face
    double v1[3], v2[3], v3[3];
    
    // Vertex 1
    v1[0] = the_face.x1;
    v1[1] = the_face.y1;
    v1[2] = the_face.z1;

    // Vertex 2
    v2[0] = the_face.x2;
    v2[1] = the_face.y2;
    v2[2] = the_face.z2;

    // Vertex 3
    v3[0] = the_face.x3;
    v3[1] = the_face.y3;
    v3[2] = the_face.z3;

    double N[3], D;
    calc_plane_coefficients(v1,v2,v3,N,D);

    double seg_rq[3];
    calc_subtract_vector(x_new,x_prox,seg_rq);

    double num, den, t;
    num = D - calc_dot_product(x_prox,N);
    den = calc_dot_product(seg_rq,N);

    // Case 1: Segment is parallel to the plane
    if (den == 0.0)
    {
        // Within the plane
        if (num == 0.0)
            return true;
        // Outside the plane
        else
            return false;
    }
    else
    {
        t = num / den;
    }

    // Case 2: Intersection point is between the endpoints
    if ( (t > 0.0) && (t < 1.0) )
        return true;
    // Case 3: Intersection point is the endpoint 'q'
    else if (t == 0.0)
        return true;
    // Case 4: Intersection point is the endpoint 'r'
    else if (t == 1.0)
        return true;
    // Case 5: There is NO intersection
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

void write_to_vtk_iteration (struct cco_network *the_network)
{
    uint32_t num_points = the_network->point_list->num_nodes;
    uint32_t num_terminals = the_network->num_terminals;
    struct point_list *p_list = the_network->point_list;
    struct point_node *p_tmp = p_list->list_nodes;
    uint32_t num_segments = the_network->segment_list->num_nodes;
    struct segment_list *s_list = the_network->segment_list;
    struct segment_node *s_tmp = s_list->list_nodes;

    char filename[MAX_FILENAME_SIZE];
    sprintf(filename,"output/cco_tree_cm_iter_%d.vtk",num_terminals);
    FILE *file = fopen(filename,"w+");

    // Write the header
    fprintf(file,"# vtk DataFile Version 3.0\n");
    fprintf(file,"Tree\n");
    fprintf(file,"ASCII\n");
    fprintf(file,"DATASET POLYDATA\n");
    fprintf(file,"POINTS %u float\n",num_points);
    while (p_tmp != NULL)
    {
        fprintf(file,"%g %g %g\n",p_tmp->value->x,p_tmp->value->y,p_tmp->value->z);
        p_tmp = p_tmp->next;
    }

    fprintf(file,"LINES %u %u\n",num_segments,num_segments*3);
    while (s_tmp != NULL)
    {
        fprintf(file,"2 %u %u\n",s_tmp->value->src->id,s_tmp->value->dest->id);
        s_tmp = s_tmp->next;
    }
    fprintf(file,"CELL_DATA %u\n",num_segments);
    fprintf(file,"SCALARS radius float\n");
    fprintf(file,"LOOKUP_TABLE default\n");
    s_tmp = s_list->list_nodes;
    while (s_tmp != NULL)
    {
        fprintf(file,"%g\n",s_tmp->value->radius * 1000.0);
        s_tmp = s_tmp->next;
    }
    fclose(file);
}
