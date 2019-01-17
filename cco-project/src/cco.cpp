#include "../include/cco.h"

Point::Point () { }

Point::Point (const double x, const double y, const double z)
{
    this->x = x;
    this->y = y;
    this->z = z;
}

Segment::Segment () { }

Segment::Segment (Point *p1, Point *p2,\
                int index_source, int index_destination,\
                Segment *left, Segment *right, Segment *parent,\
                const double Q, const double p)
{
    this->p1 = p1;
    this->p2 = p2;
    this->length = calc_size_segment(p1,p2);
    this->radius = calc_poisseulle(Q,p,this->length);
    this->left = left;
    this->index_source = index_source;
    this->index_destination = index_destination;
    this->right = right;
    this->parent = parent;
    this->type = get_segment_type();
    calc_bifurcation_ratio();
    this->Q = Q;
    this->p = p;
}

double Segment::calc_dproj (const Point p)
{
    Point *distal = this->p1;
    Point *proximal = this->p2;

    double dot_product = ((proximal->x - distal->x)*(p.x - distal->x)) +\
                         ((proximal->y - distal->y)*(p.y - distal->y)) +\
                         ((proximal->z - distal->z)*(p.z - distal->z));

    return dot_product*pow(this->length,-2.0);
}

double Segment::calc_dortho (const Point p)
{
    Point *distal = this->p1;
    Point *proximal = this->p2;

    double dot_product = fabs(((-proximal->y + distal->y)*(p.x - distal->x)) +\
                         ((proximal->x - distal->x)*(p.y - distal->y)) +\
                         ((proximal->z - distal->z)*(p.z - distal->z)));

    return dot_product*pow(this->length,-1.0);
}

double Segment::calc_dend (const Point p)
{
    Point *distal = this->p1;
    Point *proximal = this->p2;

    double d_distal = pow( pow(p.x-distal->x,2) + pow(p.y-distal->y,2) + pow(p.z-distal->z,2) ,0.5);
    double d_proximal = pow( pow(p.x-proximal->x,2) + pow(p.y-proximal->y,2) + pow(p.z-proximal->z,2) ,0.5);

    return min(d_distal,d_proximal);
}

int Segment::get_segment_type ()
{
    // Root
    if (this->parent == NULL)
        return ROOT_STICKER;
    // Terminal
    else if (this->left == NULL && this->right == NULL)
        return TERMINAL_STICKER;
    // Bifurcation
    else
        return BIFURCATION_STICKER;
}

// Calculating the bifurcation ratio following equation (1) from the paper
void Segment::calc_bifurcation_ratio ()
{
    if (this->type == TERMINAL_STICKER)
    {
        this->beta_l = UNDEFINED_VALUE;
        this->beta_r = UNDEFINED_VALUE;
    }
    else
    {
        if (this->left != NULL)
            this->beta_l = this->left->radius / this->radius;
        else
            this->beta_l = UNDEFINED_VALUE;

        if (this->right != NULL)
            this->beta_r = this->right->radius / this->radius;
        else
            this->beta_r = UNDEFINED_VALUE;
    }
}

CCO_Network::CCO_Network (User_Options *options)
{
    num_terminals = 0;
    Q_perf = options->Q_perf;
    p_perf = options->p_perf;
    r_perf = options->r_perf;
    N_term = options->N_term;
}

// Area of the supporting circle
double CCO_Network::calc_dthreashold (const double radius, const int num_term)
{
    return sqrt(M_PI*radius*radius/num_term);
}

void CCO_Network::grow_tree ()
{
    make_root();

    // Main iteration loop
    while (num_terminals <= this->N_term)
    {
        generate_new_terminal();

        num_terminals++;
    }
}

void CCO_Network::make_root ()
{
    
    // Calculating the radius of the first microcirculatory black-box (Nterm = 1 -> root)
    double r_supp = sqrt(Q_perf / M_PI);

    Point proximal(0,0,0);
    Point distal;
    
    // Calculate the distal position of the root inside the circle with radius r_supp
    srand(time(NULL));
    generate_point_inside_circle(&distal,r_supp);
    
    // Insert the points into the array of Point
    points.push_back(proximal);
    points.push_back(distal);
    print_points();

    // Build and insert the root Segment into the array of Segments
    Segment root(&points[0],&points[1],0,1,NULL,NULL,NULL,Q_perf,p_perf);
    segments.push_back(root);
    print_segments();

    num_terminals++;
}

void CCO_Network::calc_middle_segment (Point *p, const Segment segment)
{
    p->x = (segment.p1->x + segment.p2->x) / 2.0;
    p->y = (segment.p1->y + segment.p2->y) / 2.0; 
}

bool CCO_Network::has_collision (const Point p, const unsigned iconn_index)
{
    // Calculate the middle point of each segment
    Point middle_point;

    for (unsigned int i = 0; i < segments.size(); i++)
    {
        if (i != iconn_index)
        {
            Point *p1 = segments[i].p1;
            Point *p2 = segments[i].p2;
            calc_middle_segment(&middle_point,segments[i]);

            double denominator = ((p2->x - p1->x) * (p.y - middle_point.y)) - ((p2->y - p1->y) * (p.x - middle_point.x));
            double numerator1 = ((p1->y - middle_point.y) * (p.x - middle_point.x)) - ((p1->x - middle_point.x) * (p.y - middle_point.y));
            double numerator2 = ((p1->y - middle_point.y) * (p2->x - p1->x)) - ((p1->x - middle_point.x) * (p2->y - p1->y));
            
            // Detect coincident lines (has a problem, read below)
            if (denominator == 0) return numerator1 == 0 && numerator2 == 0;

            double r = numerator1 / denominator;
            double s = numerator2 / denominator;

            bool intersect = (r >= 0 && r <= 1) && (s >= 0 && s <= 1);

            if (intersect)
            {
                printf("[-] Intersection between segments!\n");
                return true;
            }
        }
    }
    
    return false;
}

void CCO_Network::generate_new_terminal ()
{
    Point new_point;
    bool point_is_ok;
    int toss;

    // Calculating the radius of the black-box
    double r_supp = sqrt(Q_perf / (num_terminals*M_PI));

    // Find the best segment to make the connection
    for (unsigned int i = 0; i < segments.size(); i++)
    {
        point_is_ok = false;
        toss = 0;

        // Calculate threashold distance
        double d_threash = calc_dthreashold(r_supp,num_terminals);

        while (!point_is_ok)
        {
            // Generate the position where the new terminal will be
            generate_point_inside_circle(&new_point,r_supp);

            // Calculate the projections
            double d_proj = segments[i].calc_dproj(new_point);

            double d_crit;
            if (d_proj >= 0 && d_proj <= 1)
                d_crit = segments[i].calc_dortho(new_point);
            else
                d_crit = segments[i].calc_dend(new_point);
            
            if (d_crit > d_threash && !has_collision(new_point,i))
            {
                printf("[+] Making new segment! Number of tosses = %d || d_threash = %.2lf\n",toss,d_threash);
                printf("[New Point] ");
                print_point(new_point);
                point_is_ok = true;
            }
            else
            {
                if (toss > N_toss)
                {
                    printf("[!] Reducing dthreash! Before = %.2lf || Now = %.2lf \n",d_threash,d_threash*0.9);
                    d_threash *= 0.9;
                    toss = 0;
                }
                //printf("[Toss %d] Making a new toss\n",toss);
                toss++;
            }
        }
        
        // Insert the new point into the array of points
        points.push_back(new_point);


    }

}

bool is_inside_circle (Point *p, const Point c, const double radius)
{
    double d = sqrt(pow(p->x - c.x,2) + pow(p->y - c.y,2) + pow(p->z - c.z,2));

    return (d <= radius) ? true : false;
}

void generate_point_inside_circle (Point *p, const double radius)
{
    // Center of the perfusion circle
    Point center; 
    center.x = 0.0;
    center.y = -radius;
    center.z = 0.0;

    double rand_number;
    do
    {
        rand_number = (double)rand() / (double)RAND_MAX;
        p->x = -radius + 2.0*rand_number*radius;

        rand_number = (double)rand() / (double)RAND_MAX;
        p->y = -2.0*rand_number*radius;

        p->z = 0.0;
         
    }while (!is_inside_circle(p,center,radius));
}

void CCO_Network::print_segments ()
{
    printf("SEGMENTS\n");
    printf("%s\n",PRINT_LINE);
    for (unsigned int i = 0; i < segments.size(); i++)
    {
        printf("%d - ",i);
        printf("(%d,%d) -- Length = %.2lf -- Radius = %.2lf -- Beta_l = %.2lf -- Beta_r = %.2lf\n",\
                                            segments[i].index_source,\
                                            segments[i].index_destination,\
                                            segments[i].length,segments[i].radius,\
                                            segments[i].beta_l,segments[i].beta_r);
    }
    printf("%s\n",PRINT_LINE);
}

void CCO_Network::print_points ()
{
    printf("POINTS\n");
    printf("%s\n",PRINT_LINE);
    for (unsigned int i = 0; i < points.size(); i++)
    {
        printf("%d - ",i);
        print_point(points[i]);
    }
    printf("%s\n",PRINT_LINE);
}

void CCO_Network::write_to_vtk ()
{
    FILE *file = fopen("cco_tree.vtk","w+");

    fprintf(file,"# vtk DataFile Version 3.0\n");
    fprintf(file,"Tree\n");
    fprintf(file,"ASCII\n");
    fprintf(file,"DATASET POLYDATA\n");
    fprintf(file,"POINTS %lu float\n",points.size());
    for (unsigned int i = 0; i < points.size(); i++)
        fprintf(file,"%g %g %g\n",points[i].x,points[i].y,points[i].z);
    fprintf(file,"LINES %lu %lu\n",segments.size(),segments.size()*3);
    for (unsigned int i = 0; i < segments.size(); i++)
        fprintf(file,"2 %d %d\n",segments[i].index_source,segments[i].index_destination);

    fclose(file);
}

void print_point (const Point p)
{
    printf("(%.2lf,%.2lf,%.2lf)\n",p.x,p.y,p.z);
}

double calc_size_segment (const Point *p1, const Point *p2)
{
    return sqrt(pow(p1->x-p2->x,2) + pow(p1->y-p2->y,2) + pow(p1->z-p2->z,2));
}

// Poisseulle's Law:
// Returns the radius for a given flux Q, pressure p and length l
double calc_poisseulle (const double Q, const double p, const double l)
{
    return pow( (8.0*Q*l*ETA)/(p*M_PI) , 0.25 );
}