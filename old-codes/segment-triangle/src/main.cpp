// Author: Lucas Berg
// Program that checks collision between a given segment and triangle.

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

// Calculate the volume of a tethahedron
double calculate_volume_sign (Point *q, Point *a, Point *b, Point *r)
{
    double vol;
    double A[3], B[3], C[3];

    A[0] = q->x - r->x;
    A[1] = q->y - r->y;
    A[2] = q->z - r->z;
    B[0] = a->x - r->x;
    B[1] = a->y - r->y;
    B[2] = a->z - r->z;
    C[0] = b->x - r->x;
    C[1] = b->y - r->y;
    C[2] = b->z - r->z;

    vol = A[0] * (B[1]*C[2] - B[2]*C[1]) +\
          A[1] * (B[2]*C[0] - B[0]*C[2]) +\
          A[2] * (B[0]*C[1] - B[1]*C[0]);

    if (vol > 0.5)          return 1;
    else if (vol < -0.5)    return -1;
    else                    return 0;
}


char check_segment_triangle_intersection (vector<Segment> triangle, vector<Segment> segment)
{
    double vol[3];

    Point *a = triangle[0].p1;
    Point *b = triangle[1].p1;
    Point *c = triangle[2].p1; 
    Point *q = segment[0].p1;
    Point *r = segment[0].p2;

    vol[0] = calculate_volume_sign(q,a,b,r);
    vol[1] = calculate_volume_sign(q,b,c,r);
    vol[2] = calculate_volume_sign(q,c,a,r);

    // Same sign: Segment intersects interior of triangle
    if ( (vol[0] > 0) && (vol[1] > 0) && (vol[2] > 0) ||\
         (vol[0] < 0) && (vol[1] < 0) && (vol[2] < 0) )
        return 'f';

    // Opposite sign: no intersection between segment and triangle
    if ( (vol[0] > 0) || (vol[1] > 0) || (vol[2] > 0) &&\
         (vol[0] < 0) || (vol[1] < 0) || (vol[2] < 0) ) 
        return '0';
    
    else if ( (vol[0] == 0.0) && (vol[1] == 0.0) && (vol[2] == 0.0))
    {
        fprintf(stderr,"[-] ERROR!\n");
        exit(EXIT_FAILURE);
    }

    // Two zeros: Segment intersects vertex
    else if ( (vol[0] == 0.0) && (vol[1] == 0.0) ||\
              (vol[0] == 0.0) && (vol[2] == 0.0) ||\
              (vol[1] == 0.0) && (vol[2] == 0.0))
        return 'v';

    // One zero: Segment intersects edge
    else if ( (vol[0] == 0.0) || (vol[1] == 0.0) || (vol[2] == 0.0) )
        return 'e';
    
    else
    {
        fprintf(stderr,"[-] ERROR!\n");
        exit(EXIT_FAILURE);
    } 

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

    char ret = check_segment_triangle_intersection(triangle,segment);
    switch (ret)
    {
        case 'v': printf("[+] The segment includes a vertex of T -- code 'v'\n");
                  break;
        case 'e': printf("[+] The segment includes a point in the interior of an edge of T -- code 'e'\n");
                  break;
        case 'f': printf("[+] The segment includes a point in the interior of a face of T -- code 'f'\n");
                  break;
        case '0': printf("[-] Segment does not intersect the triangle T -- code '0'\n");
                  break;
    }

    write_to_vtk(points,triangle,segment);

    return 0;
}