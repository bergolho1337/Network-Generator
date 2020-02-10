// Author: Lucas Berg

#include <iostream>
#include <cmath>
#include <cstring>
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

class Face
{
public:
    Point *v1;
    Point *v2;
    Point *v3;
    double normal[3];
public:
    Face (Point *a, Point *b, Point *c, double n[])
    {
        this->v1 = a;
        this->v2 = b;
        this->v3 = c;
        for (uint32_t i = 0; i < 3; i++)
            this->normal[i] = n[i];
    };
    void print ()
    {
        printf("\tVertex = [%lf %lf %lf] [%lf %lf %lf] [%lf %lf %lf]\n",this->v1->x,this->v1->y,this->v1->z,\
                                                        this->v2->x,this->v2->y,this->v2->z,\
                                                        this->v3->x,this->v3->y,this->v3->z);
        printf("\tNormal = (%lf %lf %lf)\n",this->normal[0],this->normal[1],this->normal[2]);
    };
};

void read_face (FILE *file, vector<Face> &faces)
{
    char str[200];
    double n[3];
    double a[3], b[3], c[3];

    // Read normal vector
    fscanf(file,"%s",str);
    fscanf(file,"%lf %lf %lf",&n[0],&n[1],&n[2]);

    // Read vertex
    fscanf(file,"%s",str); fscanf(file,"%s",str);
    fscanf(file,"%s",str);
    fscanf(file,"%lf %lf %lf",&a[0],&a[1],&a[2]);
    //printf("%g %g %g\n",a[0],a[1],a[2]);
    fscanf(file,"%s",str);
    fscanf(file,"%lf %lf %lf",&b[0],&b[1],&b[2]);
    //printf("%g %g %g\n",b[0],b[1],b[2]);
    fscanf(file,"%s",str);
    fscanf(file,"%lf %lf %lf",&c[0],&c[1],&c[2]);
    //printf("%g %g %g\n\n",c[0],c[1],c[2]);
    fscanf(file,"%s",str); fscanf(file,"%s",str);

    Point *v1 = new Point(0,a[0],a[1],a[2]);
    Point *v2 = new Point(1,b[0],b[1],b[2]);
    Point *v3 = new Point(2,c[0],c[1],c[2]);
    Face new_face(v1,v2,v3,n);

    faces.push_back(new_face);
}

void read_faces (const char *filename, vector<Face> &faces)
{
    char str[200];
    FILE *file = fopen(filename,"r");

    while (fscanf(file,"%s",str) != EOF)
    {
        if (strcmp(str,"facet") == 0)
            read_face(file,faces);
    }
    fclose(file);
}

void print_points (vector<Point> points)
{
    for (int i = 0; i < (int)points.size(); i++)
        points[i].print();
}

void print_faces (vector<Face> faces)
{
    for (int i = 0; i < (int)faces.size(); i++)
    {
        printf("[Face %d]\n",i);
        faces[i].print();
    }
}

int main (int argc, char *argv[])
{
    if (argc-1 != 1)
    {
        printf("Usage:> %s <input_filename>\n",argv[0]);
        exit(EXIT_FAILURE);   
    }

    vector<Face> faces;
    read_faces(argv[1],faces);
    print_faces(faces);

    return 0;
}