#include "../include/cco.h"

double D_CRIT;

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
    return sqrt( (double)(M_PI*radius*radius) / (double)(num_term) );
}

double CCO_Network::calc_dproj (const double pos[], const int iconn_index)
{
    int src = segments[iconn_index].src;
    int dest = segments[iconn_index].dest;
    double length = segments[iconn_index].length;

    Point distal = points[dest];
    Point proximal = points[src];

    double dot_product = ((proximal.x - distal.x)*(pos[0] - distal.x)) +\
                         ((proximal.y - distal.y)*(pos[1] - distal.y)) +\
                         ((proximal.z - distal.z)*(pos[2] - distal.z));

    return dot_product*pow(length,-2.0);
}

double CCO_Network::calc_dortho (const double pos[], const int iconn_index)
{
    int src = segments[iconn_index].src;
    int dest = segments[iconn_index].dest;
    double length = segments[iconn_index].length;

    Point distal = points[dest];
    Point proximal = points[src];

    double dot_product = ((-proximal.y + distal.y)*(pos[0] - distal.x)) +\
                         ((proximal.x - distal.x)*(pos[1] - distal.y)) +\
                         ((proximal.z - distal.z)*(pos[2] - distal.z));

    return fabs(dot_product)*pow(length,-1.0);
}

double CCO_Network::calc_dend (const double pos[], const int iconn_index)
{
    int src = segments[iconn_index].src;
    int dest = segments[iconn_index].dest;

    Point distal = points[dest];
    Point proximal = points[src];

    double d_distal = sqrt( pow(pos[0]-distal.x,2) + pow(pos[1]-distal.y,2) + pow(pos[2]-distal.z,2));
    double d_proximal = sqrt( pow(pos[0]-proximal.x,2) + pow(pos[1]-proximal.y,2) + pow(pos[2]-proximal.z,2));

    return min(d_distal,d_proximal);    
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

bool CCO_Network::is_inside_perfusion_area (const double pos[], const double radius)
{
    double d = sqrt(pow(pos[0] - center.x,2) + pow(pos[1] - center.y,2) + pow(pos[2] - center.z,2));

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

void CCO_Network::generate_point_inside_perfusion_area (double pos[], const double radius)
{
    double rand_number;
    do
    {
        rand_number = (double)rand() / (double)RAND_MAX;
        pos[0] = -radius + 2.0*rand_number*radius;

        rand_number = (double)rand() / (double)RAND_MAX;
        pos[1] = -2.0*rand_number*radius;

        pos[2] = 0.0;
         
    }while (!is_inside_perfusion_area(pos,radius));
}

void CCO_Network::make_root ()
{

    // Calculating the radius of the first microcirculatory black-box (Nterm = 1 -> root)
    this->r_supp = sqrt(Q_perf / (M_PI));

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
    srand(300);
    //srand(time(NULL));
    generate_point_inside_perfusion_area(&distal,r_supp);
    
    // Insert the point into the array of Points
    points.push_back(distal);
    //print_points();

    // Build and insert the root Segment into the array of Segments
    Segment root(&points[0],&points[1],NIL,NIL,NIL,Q_perf,p_perf);
    segments.push_back(root);
    print_segments();

    // TODO
    // calculate radius of the root

    num_terminals++;

    // TODO
    // Incremetar o Ndist no segmento apos a insercao

}

void CCO_Network::read_cloud_points (const char filename[], vector<Point> &cloud_points)
{
    FILE *file = fopen(filename,"r");

    unsigned int num_points;
    fscanf(file,"%u",&num_points);

    for (unsigned int i = 0; i < num_points; i++)
    {
        double pos[3];
        fscanf(file,"%lf %lf %lf",&pos[0],&pos[1],&pos[2]);

        Point p(i,pos[0],pos[1],pos[2]);
        cloud_points.push_back(p);
    }

    fclose(file);
}

void CCO_Network::grow_tree ()
{
    //test1();
    //test2();
    //test3();
    //test4();
    //test5();
    test6();

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

bool CCO_Network::collision_detect (Point p1, Point p2, Point p3, Point p4)
{
    double denominator = ((p2.x - p1.x) * (p3.y - p4.y)) - ((p2.y - p1.y) * (p3.x - p4.x));
    double numerator1 = ((p1.y - p4.y) * (p3.x - p4.x)) - ((p1.x - p4.x) * (p3.y - p4.y));
    double numerator2 = ((p1.y - p4.y) * (p2.x - p1.x)) - ((p1.x - p4.x) * (p2.y - p1.y));
    
    // Detect coincident lines (has a problem, read below)
    if (denominator == 0) return numerator1 == 0 && numerator2 == 0;

    double r = numerator1 / denominator;
    double s = numerator2 / denominator;

    return (r >= 0 && r <= 1) && (s >= 0 && s <= 1);

}

bool CCO_Network::has_collision (const double pos[], const int iconn)
{
    // Create the new point
    Point new_point(0,pos[0],pos[1],pos[2]);

    // Calculate the middle point of the connection segment
    double middle_pos[3];
    calc_middle_segment(middle_pos,&segments[iconn]);
    Point middle_point(1,middle_pos[0],middle_pos[1],middle_pos[2]);
  
    printf("[+] Trying connection with segment %u\n",iconn);
    // Then, check for collision with all segments of the tree, except iconn,
    for (int i = 0; i < (int)segments.size(); i++)
    {
        if (i != iconn)
        {
            printf("\t[!] Checking collison between segment %d\n",i);

            // Get the reference to the points of the current segment
            Point src = points[segments[i].src];
            Point dest = points[segments[i].dest];

            bool intersect = collision_detect(middle_point,new_point,src,dest);
            if (intersect)
            {
                printf("\t[-] ERROR! Intersection with segment %d !\n",i);
                return true;
            }
        }
    }
    return false;
}

int CCO_Network::connection_search (const double pos[])
{
    for (int i = 0; i < (int)segments.size(); i++)
    {
        if (!has_collision(pos,i))
        {
            printf("[!] Making connection with segment %u\n",i);
            return i;
        }
    }
    printf("[!] Collision with all segments!\n");
    return -1;
}

int CCO_Network::connection_search_closest (const double pos[])
{
    int min_index = -1;
    double min_d = DBL_MAX;

    for (int i = 0; i < (int)segments.size(); i++)
    {
        if (!has_collision(pos,i))
        {
            double new_pos[3];
            calc_middle_segment(new_pos,&segments[i]);

            double d = calc_euclidean_dist(pos,new_pos);
            if (d < min_d)
            {
                min_d = d;
                min_index = i;
            }
        }
    }
    if (min_index == -1)
        printf("[!] Collision with all segments!\n");
    return min_index;
}

bool CCO_Network::distance_criterion (const double pos[], const int iconn_index, const double d_threash)
{
    double d_proj = calc_dproj(pos,iconn_index);
    
    double d_crit;
    if (d_proj >= 0 && d_proj <= 1)
        d_crit = calc_dortho(pos,iconn_index);
    else
        d_crit = calc_dend(pos,iconn_index);
    
    if (d_crit > d_threash)
        return true;
    else
        return false;

}

int CCO_Network::connection_search_paper (const double pos[], const double d_threash)
{

    for (int i = 0; i < (int)segments.size(); i++)
    {
        // The point does NOT attend the distance criterion
        if (!distance_criterion(pos,i,d_threash))
        {
            return -1;    
        }
    }
    // The new point has passed the distance criterion for all the segments
    return 1;
}


// OK
void CCO_Network::build_segment (double new_pos[])
{
    unsigned int iconn_index = connection_search(new_pos);
    Segment *iconn = &segments[iconn_index];

    double pos[3];
    unsigned int middle_point_index = points.size(); 
    calc_middle_segment(pos,iconn);
    Point middle_point(middle_point_index,pos[0],pos[1],pos[2]);
    points.push_back(middle_point);

    // Save iconn data
    //unsigned int iconn_parent = iconn->parent;
    //unsigned int iconn_left = iconn->left;
    //unsigned int iconn_right = iconn->right;
    //unsigned int iconn_src = iconn->src;
    //unsigned int iconn_dest = iconn->dest;
    Segment *iconn_left = NULL;
    Segment *iconn_right = NULL;
    if (iconn->left != NIL)
        iconn_left = &segments[iconn->left];
    if (iconn->right != NIL)
        iconn_right = &segments[iconn->right];

    // Create ibiff
    unsigned int ibiff_index = segments.size();
    Segment ibiff(&points[middle_point_index],&points[iconn->dest],\
                iconn->left,iconn->right,iconn_index,\
                Q_perf,p_perf);
    segments[iconn_index].left = ibiff_index;
    segments[iconn_index].dest = middle_point_index;

    if (iconn_left != NULL)
        iconn_left->parent = ibiff_index;

    segments.push_back(ibiff);

    // TODO: Usar lista encadeada
    Segment *ibiff_aux = &segments[ibiff_index];
    ibiff_aux->ndist = segments[iconn_index].ndist;

    // Create inew
    unsigned int new_point_index = points.size();
    Point new_point(new_point_index,new_pos[0],new_pos[1],new_pos[2]);
    points.push_back(new_point);
    unsigned int inew_index = segments.size();
    Segment inew(&points[middle_point_index],&points[new_point_index],\
                NIL,NIL,iconn_index,\
                Q_perf,p_perf);
    segments[iconn_index].right = inew_index;

    segments[iconn_index].ndist = segments[ibiff_index].ndist + 1;
    
    if (iconn_right != NULL)
        iconn_right->parent = ibiff_index;

    segments.push_back(inew);

    // Update ndist until root
    int curr_index = segments[iconn_index].parent;
    Segment *icurr = &segments[curr_index];
    while (curr_index != NIL)
    {
        int left_index = segments[curr_index].left;
        int right_index = segments[curr_index].right;

        icurr->ndist = segments[left_index].ndist + segments[right_index].ndist;

        curr_index = segments[curr_index].parent;
        icurr = &segments[curr_index];
    }

    printf("---> Point (%lf,%lf,%lf)\n",new_pos[0],new_pos[1],new_pos[2]);
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

void CCO_Network::build_segment (const unsigned int j, double new_pos[])
{
    Segment *iconn = &segments[j];

    double pos[3];
    unsigned int middle_point_index = points.size(); 
    calc_middle_segment(pos,iconn);
    Point middle_point(middle_point_index,pos[0],pos[1],pos[2]);
    points.push_back(middle_point);

    // Save iconn data
    //unsigned int iconn_parent = iconn->parent;
    //unsigned int iconn_left = iconn->left;
    //unsigned int iconn_right = iconn->right;
    //unsigned int iconn_src = iconn->src;
    //unsigned int iconn_dest = iconn->dest;

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

    //print_segments();

}

// Update the numbers of the points
void CCO_Network::update_points (const unsigned int index)
{
    for (unsigned int i = 0; i < points.size(); i++)
    {
        if (points[i].id >= index)
            points[i].id--;
    }
}

// Update the numbers of the segment structure
void CCO_Network::update_segments (const int s_index, const unsigned int p_index)
{
    for (unsigned int i = 0; i < segments.size(); i++)
    {
        // Point indexes
        if (segments[i].src >= p_index)
            segments[i].src--;
        if (segments[i].dest >= p_index)
            segments[i].dest--;

        // Segment indexes
        if (segments[i].left >= s_index)
            segments[i].left--;
        if (segments[i].right >= s_index)
            segments[i].right--;
    }
}

void CCO_Network::destroy_offspring (const int s_index)
{
    // Get the reference to the parent Segment
    int parent_index = segments[s_index].parent;

    // Find out if the offspring to be destroy is on the 'left' or 'right'
    if (segments[parent_index].left == s_index)
        segments[parent_index].left = NIL;
    else
        segments[parent_index].right = NIL;
}

void CCO_Network::destroy_segment (const int iconn_index)
{
    
    unsigned int distal_point_index = segments[iconn_index].dest;
    points.erase(points.begin() + distal_point_index);
    update_points(distal_point_index);

    destroy_offspring(iconn_index);
    segments.erase(segments.begin() + iconn_index);
    update_segments(iconn_index,distal_point_index); 
}

void CCO_Network::generate_new_terminal_old ()
{
    // Generate the terminal position inside the perfusion area
    double pos[3];
    generate_point_inside_perfusion_area(pos,r_supp);

    // Conection search
    //int iconn_index = connection_search(pos);
    int iconn_index = connection_search_closest(pos);
    build_segment(iconn_index,pos);
}

int CCO_Network::check_collisions (const double pos[])
{
    for (int i = 0; i < (int)segments.size(); i++)
    {
        if (has_collision(pos,i))
            return -1;
    }
    return 0;
}

int CCO_Network::find_most_distant_segment (const double pos[])
{
    double max_dist = DBL_MIN;
    int max_dist_index = -1;

    for (int i = 0; i < (int)segments.size(); i++)
    {
        double d_proj = calc_dproj(pos,i);
    
        double d_crit;
        if (d_proj >= 0 && d_proj <= 1)
            d_crit = calc_dortho(pos,i);
        else
            d_crit = calc_dend(pos,i);
        

        if (d_crit > max_dist)
        {
            max_dist = d_crit;
            max_dist_index = i;
        }
    }
    return max_dist_index;
}

void CCO_Network::generate_new_terminal ()
{
    int iconn_index;
    int ret;
    bool point_is_ok = false;
    int tosses = 0;
    double pos[3];
    double d_threash = calc_dthreashold(r_supp,num_terminals);

    while (!point_is_ok)
    {
        // Generate the terminal position inside the perfusion area
        generate_point_inside_perfusion_area(pos,r_supp);

        // Check the distance criterion for this point
        ret = connection_search_paper(pos,d_threash);

        // If the point has pass the distance criterion for all the segments then 
        // we need to check collisions with other segments
        if (ret != NIL)
        {
            ret = check_collisions(pos);
            // Distance criterion is ok and there are no collisions
            if (ret != NIL)
                point_is_ok = true;
        }
        // If the point does not attend the distance criterion or if there is a collision
        // we need to choose another point.
        else
        {
            tosses++;
            if (tosses > N_toss)
            {
                printf("[!] Reducing dthreash! Before = %.2lf || Now = %.2lf \n",\
                        d_threash,d_threash*0.9);
                d_threash *= 0.9;
                tosses = 0;
            }
        }
    }
    // Now, the new point attend the distance criterion and 
    // does not collide with any segment of the tree

    // Find the closest segment to make the connection
    iconn_index = find_most_distant_segment(pos); 
    build_segment(iconn_index,pos);

}

void CCO_Network::generate_new_terminal_using_file (vector<Point> &cloud_points)
{
    int iconn_index;
    int ret;
    bool point_is_ok = false;
    int tosses = 0;
    double pos[3];
    double d_threash = calc_dthreashold(r_supp,num_terminals);

    while (!point_is_ok)
    {
        // Get a point from the cloud of points
        int index = rand() % cloud_points.size();
        pos[0] = cloud_points[index].x;
        pos[1] = cloud_points[index].y;
        pos[2] = cloud_points[index].z;

        // Check the distance criterion for this point
        ret = connection_search_paper(pos,d_threash);

        // If the point has pass the distance criterion for all the segments then 
        // we need to check collisions with other segments
        if (ret != NIL)
        {
            ret = check_collisions(pos);
            // Distance criterion is ok and there are no collisions
            if (ret != NIL)
            {
                point_is_ok = true;
                cloud_points.erase(cloud_points.begin()+index);
            }
                
        }
        // If the point does not attend the distance criterion or if there is a collision
        // we need to choose another point.
        else
        {
            tosses++;
            if (tosses > N_toss)
            {
                printf("[!] Reducing dthreash! Before = %.2lf || Now = %.2lf \n",\
                        d_threash,d_threash*0.9);
                d_threash *= 0.9;
                tosses = 0;
            }
        }
    }
    // Now, the new point attend the distance criterion and 
    // does not collide with any segment of the tree

    // Find the closest segment to make the connection
    iconn_index = find_most_distant_segment(pos); 
    build_segment(iconn_index,pos);

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

// Generate the CCO tree using a fix cloud of points
void CCO_Network::test6 ()
{
    make_root();

    vector<Point> cloud_points;
    read_cloud_points("cloud/cloud_points.txt",cloud_points);

    // Main iteration loop
    while (num_terminals <= this->N_term)
    {
        printf("%s\n",PRINT_LINE);
        printf("[!] Working on terminal number %d\n",num_terminals);            

        //generate_new_terminal();
        generate_new_terminal_using_file(cloud_points);

        num_terminals++;

        printf("%s\n",PRINT_LINE);
    }
}

// Generate the CCO tree by randomly sort the points inside the perfusion area
void CCO_Network::test5 ()
{
    make_root();

    // Main iteration loop
    while (num_terminals <= this->N_term)
    {
        printf("%s\n",PRINT_LINE);
        printf("[!] Working on terminal number %d\n",num_terminals);            

        generate_new_terminal();

        num_terminals++;

        printf("%s\n",PRINT_LINE);
    }
}

// This test the collision detection function
void CCO_Network::test4 ()
{
    // Root
    Point A(0,0,0,0);
    Point B(1,0,-4,0);
    points.push_back(A);
    points.push_back(B);

    Segment s1(&A,&B,NIL,NIL,NIL,Q_perf,p_perf);
    segments.push_back(s1);

    // First insert
    double pos1[3] = {-2,-4,0};
    build_segment(pos1);
    printf("Segment 1 has been inserted sucessfully!\n");

    // Second insert
    double pos2[3] = {-1,-4,0};
    build_segment(pos2);
    printf("Segment 2 has been inserted sucessfully!\n");

    // Third insert
    double pos3[3] = {2,-2,0};
    build_segment(pos3);
    printf("Segment 3 has been inserted sucessfully!\n");

    // Forth insert
    double pos4[3] = {-0.5,-5,0};
    build_segment(pos4);
    printf("Segment 4 has been inserted sucessfully!\n");

    print_segments();
}

// This test the insert and delete segment functions
void CCO_Network::test3 ()
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
    double pos2[3] = {-2,-4,0};
    build_segment(1,pos2);

    // Third segment
    double pos3[3] = {-1,-3,0};
    build_segment(2,pos3);

    destroy_segment(3);
    destroy_segment(5);

}

// This test the insert segment function
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
    double pos2[3] = {-2,-4,0};
    build_segment(1,pos2);

    // Third segment
    double pos3[3] = {-1,-3,0};
    build_segment(2,pos3);

    print_segments();

}

// This test the Point and Segment data structure
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

double calc_euclidean_dist (const double a[], const double b[])
{
    return sqrt(pow(a[0]-b[0],2) + pow(a[1]-b[1],2) + pow(a[2]-b[2],2));
}
