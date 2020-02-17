// Author: Lucas Berg
// Program that reads a STL file and store its faces and point on a vector.

#include <iostream>
#include <cmath>
#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <map>

using namespace std;

static const double dx = 12.0;
static const double dy = 12.0;

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
    bool operator <(const Point& pt) const
    {
        return (this->x < pt.x) || \
               ((!(pt.x < this->x)) && (this->y < pt.y)) || \
               ((!(pt.x < this->x)) && (!(pt.y < this->y)) && (this->z < pt.z));
    }
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

void insert_points_from_faces_to_map (vector<Face> faces, map<Point,uint32_t> &unique_points)
{
    map<Point,uint32_t>::iterator it, it2;
    uint32_t unique_num_points = 0;
    for (uint32_t i = 0; i < faces.size(); i++)
    {
        Point *v1 = faces[i].v1;
        Point *v2 = faces[i].v2;
        Point *v3 = faces[i].v3;

        it = unique_points.find(*v1);
        if (it == unique_points.end())
        {
            unique_points.insert(pair<Point,uint32_t>(*v1,unique_num_points));
            unique_num_points++;
        }
        it = unique_points.find(*v2);
        if (it == unique_points.end())
        {
            unique_points.insert(pair<Point,uint32_t>(*v2,unique_num_points));
            unique_num_points++;
        }
        it = unique_points.find(*v3);
        if (it == unique_points.end())
        {
            unique_points.insert(pair<Point,uint32_t>(*v3,unique_num_points));
            unique_num_points++;
        }
    }

    // Map points to face
    vector< vector< int > > points_to_faces;
    points_to_faces.assign(unique_points.size(),vector<int>());

    for (it = unique_points.begin(); it != unique_points.end(); ++it)
    {
        // Check if the current point matches one of the 3 points from the faces
        uint32_t cur_id = it->second;
        printf("Point %u -- (%lf %lf %lf)\n",it->second,it->first.x,it->first.y,it->first.z);

        for (uint32_t j = 0; j < faces.size(); j++)
        {
            Point *v1 = faces[j].v1;
            Point *v2 = faces[j].v2;
            Point *v3 = faces[j].v3;

            it2 = unique_points.find(*v1);
            if (it2->second == cur_id)
            {
                points_to_faces[cur_id].push_back(j);
            }
            it2 = unique_points.find(*v2);
            if (it2->second == cur_id)
            {
                points_to_faces[cur_id].push_back(j);
            }
            it2 = unique_points.find(*v3);
            if (it2->second == cur_id)
            {
                points_to_faces[cur_id].push_back(j);
            }
        }
    }
    
    for (uint32_t i = 0; i < points_to_faces.size(); i++)
    {
        printf("Point %u -- ",i);
        for (uint32_t j = 0; j < points_to_faces[i].size(); j++)
        {
            printf("%u ",points_to_faces[i][j]);
        }
        printf("\n");
    }
    //printf("Point %u -- (%lf %lf %lf)\n",it->second,it->first.x,it->first.y,it->first.z);   
} 

void draw_triangles (FILE *file, const double x, const double y)
{

    // First triangle
    fprintf(file," facet normal 0 0 1\n");
    fprintf(file,"  outer loop\n");
    fprintf(file,"   vertex %lf %lf 0\n",x,y);
    fprintf(file,"   vertex %lf %lf 0\n",x+dx,y);
    fprintf(file,"   vertex %lf %lf 0\n",x+dx,y+dy);
    fprintf(file,"  endloop\n");
    fprintf(file," endfacet\n");
    
    // Second triangle
    fprintf(file," facet normal 0 0 1\n");
    fprintf(file,"  outer loop\n");
    fprintf(file,"   vertex %lf %lf 0\n",x,y);
    fprintf(file,"   vertex %lf %lf 0\n",x,y+dy);
    fprintf(file,"   vertex %lf %lf 0\n",x+dx,y+dy);
    fprintf(file,"  endloop\n");
    fprintf(file," endfacet\n");
}

void generate_grid_stl (const uint32_t nx, const uint32_t ny)
{
    FILE *file = fopen("output.stl","w+");
    fprintf(file,"solid ascii\n");

    for (uint32_t i = 0; i < nx; i++)
    {
        double x = i*dx;
        for (uint32_t j = 0; j < ny; j++)
        {
            double y = j*dy;

            draw_triangles(file,x,y);
        }
    }

    fprintf(file,"endsolid\n");    
    fclose(file);
}

int main (int argc, char *argv[])
{
    if (argc-1 != 1)
    {
        printf("Usage:> %s <input_filename>\n",argv[0]);
        exit(EXIT_FAILURE);   
    }

    const uint32_t n = 50;
    generate_grid_stl(n,n);

    //vector<Face> faces;
    //read_faces(argv[1],faces);
    //print_faces(faces);

    //map<Point,uint32_t> unique_points;
    //insert_points_from_faces_to_map(faces,unique_points);
    
    return 0;
}