// Author: Lucas Berg
// Program that checks the collision between a segment and the plane composed by a given triangle. 

#include <iostream>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>

using namespace std;

class Point
{
public:
    int id;
    double x, y, z;
public:
    Point () {};
    Point (const int id, const double x, const double y, const double z)
    {
        this->id = id;
        this->x = x;
        this->y = y;
        this->z = z;
    };
    void print ()
    {
        printf("%d -- (%.2lf,%.2lf,%.2lf)\n",this->id,this->x,this->y,this->z);
    };
};

class Segment
{
public:
    Point *p1;
    Point *p2;
public:
    Segment (Point *a, Point *b)
    {
        this->p1 = a;
        this->p2 = b;
    };
    void print ()
    {
        printf("(%d -- %d)\n",this->p1->id,this->p2->id);
    };
};

void read_points (const char *filename, vector<Point> &points)
{
    FILE *file = fopen(filename,"r");
    double x, y, z;

    while (fscanf(file,"%lf %lf %lf",&x,&y,&z) != EOF)
    {
        Point p(points.size(),x,y,z);
        points.push_back(p);
    }
    fclose(file);
}

void print_points (vector<Point> points)
{
    for (int i = 0; i < (int)points.size(); i++)
        points[i].print();
}

void print_segments (vector<Segment> segments)
{
    for (int i = 0; i < (int)segments.size(); i++)
        segments[i].print();
}

void write_to_vtk (vector<Point> points, vector<Segment> triangle, vector<Segment> segment)
{
    uint32_t num_points = points.size();
    uint32_t num_lines = triangle.size() + segment.size();

    FILE *file = fopen("lines.vtk","w+");

    fprintf(file,"# vtk DataFile Version 3.0\n");
    fprintf(file,"Line Collision\n");
    fprintf(file,"ASCII\n");
    fprintf(file,"DATASET POLYDATA\n");
    fprintf(file,"POINTS %u float\n",num_points);
    for (unsigned int i = 0; i < num_points; i++)
        fprintf(file,"%g %g %g\n",points[i].x,points[i].y,points[i].z);
    fprintf(file,"LINES %u %u\n",num_lines,num_lines*3);
    for (unsigned int i = 0; i < triangle.size(); i++)
        fprintf(file,"2 %d %d\n",triangle[i].p1->id,triangle[i].p2->id);
    for (unsigned int i = 0; i < segment.size(); i++)
        fprintf(file,"2 %d %d\n",segment[i].p1->id,segment[i].p2->id);

    fclose(file);
}

void write_intersection_point_to_vtk (Point p)
{
    FILE *file = fopen("intersection.vtk","w+");

    fprintf(file,"# vtk DataFile Version 3.0\n");
    fprintf(file,"Intersection Point\n");
    fprintf(file,"ASCII\n");
    fprintf(file,"DATASET POLYDATA\n");
    fprintf(file,"POINTS 1 float\n");
    fprintf(file,"%g %g %g\n",p.x,p.y,p.z);
    fprintf(file,"VERTICES 1 2\n");
    fprintf(file,"1 0\n");

    fclose(file);
}

double dot_product (const double x1, const double y1, const double z1,\
                    const double x2, const double y2, const double z2)
{
    return x1*x2 + y1*y2 + z1*z2; 
}

void subtract_vector (const Point *r, const Point *q, double rq[])
{
    rq[0] = r->x - q->x;
    rq[1] = r->y - q->y;
    rq[2] = r->z - q->z;
}

void calculate_normal_vector(Point *a, Point *b, Point *c, double N[])
{
    N[0] = ( c->z - a->z ) * ( b->y - a->y ) -\
           ( b->z - a->z ) * ( c->y - a->y);
    N[1] = ( b->z - a->z ) * ( c->x - a->x ) -\
           ( b->x - a->x ) * ( c->z - a->z);
    N[2] = ( b->x - a->x ) * ( c->y - a->y ) -\
           ( b->y - a->y ) * ( c->x - a->x);
}

uint32_t calculate_plane_coefficients (vector<Segment> triangle, double N[], double &D)
{
    Point *p1 = triangle[0].p1;
    Point *p2 = triangle[1].p1;
    Point *p3 = triangle[2].p1;

    calculate_normal_vector(p1,p2,p3,N);
    D = dot_product(N[0],N[1],N[2],p1->x,p1->y,p1->z);

    // Find the largest component of N
    double largest = 0.0;
    uint32_t largest_index = 0;
    for (uint32_t i = 0; i < 3; i++)
    {
        double value = fabs(N[i]);
        if (value > largest)
        {
            largest = value;
            largest_index = i;
        }
    }

    return largest_index;
} 

char check_segment_plane_intersection (vector<Segment> triangle, vector<Segment> segment, Point &p)
{
    double segment_rq[3];   // Segment vector 
    double N[3];            // Triangle plane normal vector
    double D;               // 'D' coefficient of the plane equation
    uint32_t m;             // Index of the largest component of N
    double num, den, t;

    Point *q = segment[0].p1;
    Point *r = segment[0].p2;

    m = calculate_plane_coefficients(triangle,N,D);
    subtract_vector(r,q,segment_rq);

    num = D - dot_product(q->x,q->y,q->z,N[0],N[1],N[2]);
    den = dot_product(segment_rq[0],segment_rq[1],segment_rq[2],N[0],N[1],N[2]);

    // Case 1: Segment is parallel to the plane
    if (den == 0.0)
    {
        printf("[!] Segment is parallel to the plane\n");
        // Endpoint 'q' is on plane
        if (num == 0.0)
            return 'p';
        else
            return '0';
    }
    else
    {
        t = num / den;
    }

    // DEBUG
    //printf("N = (%g,%g,%g)\n",N[0],N[1],N[2]);
    //printf("num = %g\n",num);
    //printf("den = %g\n",den);
    //printf("t = %g\n",t);

    // Calculate the point of intersection
    p.x = q->x + t * segment_rq[0];
    p.y = q->y + t * segment_rq[1];
    p.z = q->z + t * segment_rq[2];
    printf("[+] Intersection point = (%g,%g,%g)\n",p.x,p.y,p.z);
    write_intersection_point_to_vtk(p);

    // Case 2: Intersection point is between the endpoints
    if ( (t > 0.0) && (t < 1.0) )
        return '1';
    // Case 3: Intersection point is the endpoint 'q'
    else if (t == 0.0)
        return 'q';
    // Case 4: Intersection point is the endpoint 'r'
    else if (t == 1.0)
        return 'r';
    // Case 5: There is NO intersection
    else
        return '0'; 
}

int main (int argc, char *argv[])
{
    if (argc-1 != 1)
    {
        printf("Usage:> %s <input_filename>\n",argv[0]);
        exit(EXIT_FAILURE);   
    }
    vector<Point> points;
    read_points(argv[1],points);
    print_points(points);

    vector<Segment> triangle;
    Segment s1(&points[0],&points[1]);
    triangle.push_back(s1);
    Segment s2(&points[1],&points[2]);
    triangle.push_back(s2);
    Segment s3(&points[2],&points[0]);
    triangle.push_back(s3);
    print_segments(triangle);

    vector<Segment> segment;
    Segment s4(&points[3],&points[4]);
    segment.push_back(s4);

    Point p;

    char ret = check_segment_plane_intersection(triangle,segment,p);
    switch (ret)
    {
        case 'p': printf("[+] Segment lies within the plane -- code 'p'\n");
                  break;
        case 'q': printf("[+] The (first) q endpoint is on the plane -- code 'q'\n");
                  break;
        case 'r': printf("[+] The (second) r endpoint is on the plane -- code 'r'\n");
                  break;
        case '0': printf("[-] Segment is outside the plane -- code '0'\n");
                  break;
        case '1': printf("[+] Segment intersects the plane -- code '1'\n");
                  break;
    }

    write_to_vtk(points,triangle,segment);

    return 0;
}