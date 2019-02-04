#include "../include/cco.h"

Point::Point () { }

Point::Point (const int id, const double x, const double y, const double z)
{
    this->id = id;
    this->x = x;
    this->y = y;
    this->z = z;
}

Segment::Segment () { }

Segment::Segment (Point *p1, Point *p2,\
                int left, int right, int parent,\
                const double Q, const double p)
{
    this->src = p1->id;
    this->dest = p2->id;
    this->left = left;
    this->right = right;
    this->parent = parent;
    this->Q = Q;
    this->p = p;
    this->ndist = 1;

    this->type = get_segment_type();
    this->length = calc_size_segment(p1,p2);
    this->radius = calc_poisseulle(this->Q,this->p,this->length);

    //calc_bifurcation_ratio(this);

}

// Area of the supporting circle
double CCO_Network::calc_dthreashold (const double radius, const int num_term)
{
    return sqrt(M_PI*radius*radius/num_term);
}

double Segment::calc_dproj (const Point p)
{
    /*
    Point *distal = this->p1;
    Point *proximal = this->p2;

    double dot_product = ((proximal->x - distal->x)*(p.x - distal->x)) +\
                         ((proximal->y - distal->y)*(p.y - distal->y)) +\
                         ((proximal->z - distal->z)*(p.z - distal->z));

    return dot_product*pow(this->length,-2.0);
    */
   return 0;
}

double Segment::calc_dortho (const Point p)
{
    /*
    Point *distal = this->p1;
    Point *proximal = this->p2;

    double dot_product = fabs(((-proximal->y + distal->y)*(p.x - distal->x)) +\
                         ((proximal->x - distal->x)*(p.y - distal->y)) +\
                         ((proximal->z - distal->z)*(p.z - distal->z)));

    return dot_product*pow(this->length,-1.0);
    */
   return 0;
}

double Segment::calc_dend (const Point p)
{
    /*
    Point *distal = this->p1;
    Point *proximal = this->p2;

    double d_distal = pow( pow(p.x-distal->x,2) + pow(p.y-distal->y,2) + pow(p.z-distal->z,2) ,0.5);
    double d_proximal = pow( pow(p.x-proximal->x,2) + pow(p.y-proximal->y,2) + pow(p.z-proximal->z,2) ,0.5);

    return min(d_distal,d_proximal);
    */
   return 0;
}

void Segment::add_offspring (Segment *new_segment)
{
    /*
    this->left = new_segment;
    if (this->left == NULL)
        this->left = new_segment;
    else if (this->right == NULL)
        this->right = new_segment;
    else
    {
        printf("[-] ERROR! Segment has already 2 offsprings\n");
    }
    */
} 

int Segment::get_segment_type ()
{
    // Root
    if (this->parent == NIL)
        return ROOT_STICKER;
    // Terminal
    else if (this->left == NIL && this->right == NIL)
        return TERMINAL_STICKER;
    // Bifurcation
    else
        return BIFURCATION_STICKER;
}

// Calculating the bifurcation ratio following equation (1) from the paper
void CCO_Network::calc_bifurcation_ratio (Segment *s)
{
    if (s->type == TERMINAL_STICKER)
    {
        s->beta_l = NIL;
        s->beta_r = NIL;
    }
    else
    {

        if (s->left != NIL)
        {
            Segment *left = &segments[s->left];
            s->beta_l = left->radius / s->radius;
        }
            
        else
            s->beta_l = NIL;

        if (s->right != NIL)
        {
            Segment *right = &segments[s->right];
            s->beta_r = right->radius / s->radius;
        }
            
        else
            s->beta_r = NIL;
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

bool CCO_Network::is_inside_perfusion_area (const Point *p, const double radius)
{
    double d = sqrt(pow(p->x - center.x,2) + pow(p->y - center.y,2) + pow(p->z - center.z,2));

    return (d <= radius) ? true : false;
}

void CCO_Network::generate_point_inside_perfusion_area (Point *p, const double radius)
{
    double rand_number;
    unsigned int num_points = points.size();
    do
    {
        rand_number = (double)rand() / (double)RAND_MAX;
        p->x = -radius + 2.0*rand_number*radius;

        rand_number = (double)rand() / (double)RAND_MAX;
        p->y = -2.0*rand_number*radius;

        p->z = 0.0;
         
    }while (!is_inside_perfusion_area(p,radius));
    p->id = num_points;
}

void CCO_Network::make_root ()
{

    // Calculating the radius of the first microcirculatory black-box (Nterm = 1 -> root)
    this->r_supp = sqrt(Q_perf / M_PI);

    // Set the center point of the perfusion area
    center.x = 0.0;
    center.y = -r_supp;
    center.z = 0.0;
    draw_perfusion_area(r_supp);

    // Create the proximal and distal points
    Point proximal(0,0,0,0);
    points.push_back(proximal);
    
    Point distal;
    
    // Calculate the distal position of the root inside the circle with radius r_supp
    srand(time(NULL));
    generate_point_inside_perfusion_area(&distal,r_supp);
    
    // Insert the point into the array of Points
    points.push_back(distal);
    //print_points();

    // Build and insert the root Segment into the array of Segments
    Segment root(&points[0],&points[1],NIL,NIL,NIL,Q_perf,p_perf);
    segments.push_back(root);
    print_segments();

    num_terminals++;

}

void CCO_Network::grow_tree ()
{
    //make_root();

    //test1();

    test2();

    /*
    // Main iteration loop
    while (num_terminals <= this->N_term)
    {
        printf("%s\n",PRINT_LINE);
        printf("[!] Working on terminal number %d\n",num_terminals);            

        generate_new_terminal();

        num_terminals++;
        printf("%s\n",PRINT_LINE);
    }
    */
}

void CCO_Network::draw_perfusion_area (const double radius)
{
    // Create a circle
    vtkSmartPointer<vtkRegularPolygonSource> polygonSource =
      vtkSmartPointer<vtkRegularPolygonSource>::New();
    
    polygonSource->GeneratePolygonOff(); // Outline of the circle
    polygonSource->SetNumberOfSides(50);
    polygonSource->SetRadius(radius);
    polygonSource->SetCenter(0,-radius,0);

    vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
    writer->SetFileName("perfusion_area.vtp");
    writer->SetInputConnection(polygonSource->GetOutputPort());
    writer->Write();

}




void CCO_Network::calc_middle_segment (double pos[], Segment *s)
{
    Point *p1 = &points[s->src];
    Point *p2 = &points[s->dest];
    
    pos[0] = (p1->x + p2->x) / 2.0;
    pos[1] = (p1->y + p2->y) / 2.0;
    pos[2] = (p1->z + p2->z) / 2.0;
}

// TODO: Revise
bool CCO_Network::has_collision (Point p, const unsigned iconn_index)
{
/*
    // Calculate the middle point of the connection segment
    Point middle_point;
    calc_middle_segment(&middle_point,segments[iconn_index]);
  
    // Then, check for collision with all segments of the tree, except iconn,
    for (unsigned int i = 0; i < segments.size(); i++)
    {
        if (i != iconn_index)
        {
            printf("Checking collison between segment %d -- %d\n",i,iconn_index);

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
                printf("[-] Intersection with segment %u !\n",i);
                return true;
            }
        }
    }
*/
    return false;
}

void CCO_Network::get_feasible_point (Point *p, const double radius)
{

    // Save the current number of segments
    unsigned int curr_num_segment = segments.size();

    // Calculate threashold distance
    double d_threash = calc_dthreashold(radius,num_terminals);

    bool point_is_ok = false;
    int toss = 0;

    // Until we not found a suitable point we will repeated this loop
    Point new_point;
    while (!point_is_ok)
    {
        // Generate the position where the new terminal will be
        generate_point_inside_circle(&new_point,radius);

        // Test if current point will pass the distance criterion
        for (unsigned int j = 0; j < curr_num_segment; j++)
        {
            // Calculate the projections
            double d_proj = segments[j].calc_dproj(new_point);

            double d_crit;
            if (d_proj >= 0 && d_proj <= 1)
                d_crit = segments[j].calc_dortho(new_point);
            else
                d_crit = segments[j].calc_dend(new_point);

            // TODO: Perguntar para o Rafael se o teste tem ser validado para todos os segmentos
            // The point is feasible
            if (d_crit > d_threash)
            {
                printf("[+] Making new segment! Number of tosses = %d || d_threash = %.2lf\n",toss,d_threash);
                printf("[New Point] ");
                print_point(new_point);
                point_is_ok = true;
                break;
            }
            else
            {
                if (toss > N_toss)
                {
                    printf("[!] Reducing dthreash! Before = %.2lf || Now = %.2lf \n",d_threash,d_threash*0.9);
                    d_threash *= 0.9;
                    toss = 0;
                }
                toss++;
            }            
        }

    }

    p->x = new_point.x;
    p->y = new_point.y;
    p->z = new_point.z;
}

void CCO_Network::build_segment (const unsigned int j, double new_pos[])
{
    Segment *iconn = &segments[j];

    double pos[3];
    unsigned int middle_point_index = points.size(); 
    calc_middle_segment(pos,iconn);
    Point middle_point(middle_point_index,pos[0],pos[1],pos[2]);
    points.push_back(middle_point);

    // Save iconn data
    unsigned int iconn_parent = iconn->parent;
    unsigned int iconn_left = iconn->left;
    unsigned int iconn_right = iconn->right;
    unsigned int iconn_src = iconn->src;
    unsigned int iconn_dest = iconn->dest;

    // Create ibiff
    unsigned int ibiff_index = segments.size();
    Segment ibiff(&points[middle_point_index],&points[iconn->dest],\
                iconn->left,iconn->right,j,\
                Q_perf,p_perf);
    segments[j].left = ibiff_index;
    segments[j].dest = middle_point_index;
    segments.push_back(ibiff);

    // Create inew
    unsigned int new_point_index = points.size();
    Point new_point(new_point_index,new_pos[0],new_pos[1],new_pos[2]);
    points.push_back(new_point);
    unsigned int inew_index = segments.size();
    Segment inew(&points[middle_point_index],&points[new_point_index],\
                NIL,NIL,j,\
                Q_perf,p_perf);
    segments[j].right = inew_index;
    segments.push_back(inew);

    print_segments();

}

// Construct a new segment from a Segment 'j' of the current tree
void CCO_Network::build_segment (const unsigned int j)
{
    double pos[3];

    Segment *iconn = &segments[j];
    
    calc_middle_segment(pos,iconn);
    Point middle_point(points.size(),pos[0],pos[1],pos[2]);
    points.push_back(middle_point);

    Point new_point;
    generate_point_inside_perfusion_area(&new_point,this->r_supp);
    points.push_back(new_point);
}

void CCO_Network::destroy_segment (const int j)
{

}

void CCO_Network::generate_new_terminal ()
{
/*
    // Calculating the radius of the black-box
    double radius = sqrt(Q_perf / M_PI);
    
    Point new_point;
    get_feasible_point(&new_point,radius);

    // Connection search
    unsigned int num_segments = segments.size();
    for (unsigned int j = 0; j < num_segments; j++)
    {
        if (!has_collision(new_point,j))
        {
            int ibiff_index = construct_segment(j,new_point);
        }
    }
*/
}

void CCO_Network::create_bifurcation (const int iconn_index, Point new_point)
{
    
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
        Point *p1 = &points[segments[i].src];
        Point *p2 = &points[segments[i].dest];
        printf("(%u,%u) -- Length = %.2lf -- Radius = %.2lf -- Beta_l = %.2lf -- Beta_r = %.2lf -- Q = %.2lf -- NDIST = %d\n",\
                                            p1->id,\
                                            p2->id,\
                                            segments[i].length,segments[i].radius,\
                                            segments[i].beta_l,segments[i].beta_r,\
                                            segments[i].Q,\
                                            segments[i].ndist);
        
        printf("\t");
        if (segments[i].parent != NIL)
            printf("Parent = %u",segments[i].parent);
        else
            printf("Parent = NIL");
        if (segments[i].left != NIL)
            printf(" -- Left = %u",segments[i].left);
        else
            printf(" -- Left = NIL");
        if (segments[i].right != NIL)
            printf(" -- Right = %u\n",segments[i].right);
        else
            printf(" -- Right = NIL\n");
    }
    printf("%s\n",PRINT_LINE);
}

void CCO_Network::print_points ()
{
    printf("POINTS\n");
    printf("%s\n",PRINT_LINE);
    for (unsigned int i = 0; i < points.size(); i++)
    {
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
        fprintf(file,"2 %d %d\n",segments[i].src,segments[i].dest);

    fclose(file);
}

void print_point (const Point p)
{
    printf("%u = (%.2lf,%.2lf,%.2lf)\n",p.id,p.x,p.y,p.z);
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

void CCO_Network::test2 ()
{
    // Root
    Point A(0,0,0,0);
    Point B(1,-3,-3,0);
    points.push_back(A);
    points.push_back(B);

    Segment s1(&A,&B,NIL,NIL,NIL,Q_perf,p_perf);
    segments.push_back(s1);

    // First segment
    double pos1[3] = {1,-3,0};
    build_segment(0,pos1);

    // Second segment
    //double pos2[3] = {-2,-4,0};
    //build_segment(1,pos2);

    print_segments();

}

void CCO_Network::test1 ()
{
    Point A(0,0,0,0);
    Point B(1,-3,-3,0);
    Point C(2,-1,-1,0);
    Point D(3,1,-3,0);
    Point E(4,-2,-2,0);
    Point F(5,-2,-4,0);
    Point G(6,0,-2,0);
    Point H(7,-1,-3,0);

    points.push_back(A);
    points.push_back(B);
    points.push_back(C);
    points.push_back(D);
    points.push_back(E);
    points.push_back(F);
    points.push_back(G);
    points.push_back(H);

    Segment s1(&A,&C,NIL,NIL,NIL,Q_perf,p_perf);
    Segment s2(&C,&E,NIL,NIL,NIL,Q_perf,p_perf);
    Segment s3(&C,&G,NIL,NIL,NIL,Q_perf,p_perf);
    Segment s4(&G,&H,NIL,NIL,NIL,Q_perf,p_perf);
    Segment s5(&G,&D,NIL,NIL,NIL,Q_perf,p_perf);
    Segment s6(&E,&B,NIL,NIL,NIL,Q_perf,p_perf);
    Segment s7(&E,&F,NIL,NIL,NIL,Q_perf,p_perf);

    segments.push_back(s1);
    segments.push_back(s2);
    segments.push_back(s3);
    segments.push_back(s4);
    segments.push_back(s5);
    segments.push_back(s6);
    segments.push_back(s7);

}