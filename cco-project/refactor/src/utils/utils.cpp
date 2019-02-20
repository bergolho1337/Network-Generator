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
    double r = generate_random_number();

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
    struct point *distal = s->value->src->value;
    struct point *proximal = s->value->dest->value;

    double length = euclidean_norm(proximal->x,proximal->y,proximal->z,\
                                   distal->x,distal->y,distal->z);

    double dot_product = (proximal->x - distal->x)*(pos[0] - distal->x) +\
                         (proximal->y - distal->y)*(pos[1] - distal->y);
    
    return dot_product * pow(length,-2.0);
}

double calc_dortho (struct segment_node *s, const double pos[])
{
    struct point *distal = s->value->src->value;
    struct point *proximal = s->value->dest->value;

    double length = euclidean_norm(proximal->x,proximal->y,proximal->z,\
                                   distal->x,distal->y,distal->z);

    double dot_product = (-proximal->y + distal->y)*(pos[0] - distal->x) +\
                        (proximal->x - distal->x)*(pos[1] - distal->y);

    return fabs(dot_product) * pow(length,-1.0);
}

double calc_dend (struct segment_node *s, const double pos[])
{
    struct point *distal = s->value->src->value;
    struct point *proximal = s->value->dest->value;

    double d_dist = euclidean_norm(pos[0],pos[1],pos[2],\
                                       distal->x,distal->y,distal->z);

    double d_prox = euclidean_norm(pos[0],pos[1],pos[2],\
                                    proximal->x,proximal->y,proximal->z);
    
    return std::min(d_dist,d_prox);
}

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

bool collision_detection (const double x1, const double y1, const double z1,\
                          const double x2, const double y2, const double z2,\
                          const double x3, const double y3, const double z3,\
                          const double x4, const double y4, const double z4)
{
    double denominator = ((x2 - x1) * (y3 - y4)) - ((y2 - y1) * (x3 - x4));
    double numerator1 = ((y1 - y4) * (x3 - x4)) - ((x1 - x4) * (y3 - y4));
    double numerator2 = ((y1 - y4) * (x2 - x1)) - ((x1 - x4) * (y2 - y1));
    
    // Detect coincident lines (has a problem, read below)
    if (denominator == 0) return numerator1 == 0 && numerator2 == 0;

    double r = numerator1 / denominator;
    double s = numerator2 / denominator;

    return (r >= 0 && r <= 1) && (s >= 0 && s <= 1);
}