#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>

#include <map>
#include <vector>

#define PRINT_LINE "==========================================================================================="
#define PRINT_DOTS "..........................................................................................."

class Point_Custom
{
public:
    double x, y, z;
public:
    Point_Custom () {};
    Point_Custom (const double x, const double y, const double z)
    {
        this->x = x;
        this->y = y;
        this->z = z;
    };
    void print ()
    {
        printf("(%.2lf,%.2lf,%.2lf)\n",this->x,this->y,this->z);
    };
    bool operator <(const Point_Custom& pt) const
    {
        return (this->x < pt.x) || \
               ((!(pt.x < this->x)) && (this->y < pt.y)) || \
               ((!(pt.x < this->x)) && (!(pt.y < this->y)) && (this->z < pt.z));
    }
};

class Face_Custom
{
public:
    Point_Custom *v1;
    Point_Custom *v2;
    Point_Custom *v3;
    double normal[3];
public:
    Face_Custom (Point_Custom *a, Point_Custom *b, Point_Custom *c, double n[])
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

// GLOBAL VARIABLES
std::vector<Face_Custom> mesh_faces;
std::map<Point_Custom,uint32_t> unique_points;
std::vector< std::vector< uint32_t > > points_to_faces;

void read_face (FILE *file, std::vector<Face_Custom> &faces)
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

    Point_Custom *v1 = new Point_Custom(a[0],a[1],a[2]);
    Point_Custom *v2 = new Point_Custom(b[0],b[1],b[2]);
    Point_Custom *v3 = new Point_Custom(c[0],c[1],c[2]);
    Face_Custom new_face(v1,v2,v3,n);

    faces.push_back(new_face);
}

void read_faces_from_stl (const char filename[])
{
    printf("[!] Reading STL file with the mesh nodes ...\n");

    FILE *file = fopen(filename,"r");
    if (!file)
    {
        fprintf(stderr,"[custom] ERROR ! Cannot open STL file!\n");
        exit(EXIT_FAILURE);
    }

    char str[200];
    while (fscanf(file,"%s",str) != EOF)
    {
        if (strcmp(str,"facet") == 0)
            read_face(file,mesh_faces);
    }

    fclose(file);
}

void write_points_from_faces_to_map ()
{
    printf("[!] Calculating points to faces mapping ...\n");

    std::map<Point_Custom,uint32_t>::iterator it, it2;
    uint32_t unique_num_points = 0;
    for (uint32_t i = 0; i < mesh_faces.size(); i++)
    {
        Point_Custom *v1 = mesh_faces[i].v1;
        Point_Custom *v2 = mesh_faces[i].v2;
        Point_Custom *v3 = mesh_faces[i].v3;

        it = unique_points.find(*v1);
        if (it == unique_points.end())
        {
            unique_points.insert(std::pair<Point_Custom,uint32_t>(*v1,unique_num_points));
            unique_num_points++;
        }
        it = unique_points.find(*v2);
        if (it == unique_points.end())
        {
            unique_points.insert(std::pair<Point_Custom,uint32_t>(*v2,unique_num_points));
            unique_num_points++;
        }
        it = unique_points.find(*v3);
        if (it == unique_points.end())
        {
            unique_points.insert(std::pair<Point_Custom,uint32_t>(*v3,unique_num_points));
            unique_num_points++;
        }
    }

    // Map points to face
    points_to_faces.assign(unique_points.size(),std::vector<uint32_t>());

    for (it = unique_points.begin(); it != unique_points.end(); ++it)
    {
        // Check if the current point matches one of the 3 points from the faces
        uint32_t cur_id = it->second;
        //printf("Point %u -- (%lf %lf %lf)\n",it->second,it->first.x,it->first.y,it->first.z);

        for (uint32_t j = 0; j < mesh_faces.size(); j++)
        {
            Point_Custom *v1 = mesh_faces[j].v1;
            Point_Custom *v2 = mesh_faces[j].v2;
            Point_Custom *v3 = mesh_faces[j].v3;

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

    // Write the mapping graph to a file
    // FORMAT:  The first line is the number of points
    //          Followed by the coordinates of the points
    //          Next, each line of the file represents a link between a point and face
    //          <num_points>
    //          <x> <y> <z>
    //              .
    //              .
    //              .
    //          <point_i> <face_j>
    //              .
    //              .
    //              .
    FILE *file = fopen("outputs/mapping.txt","w+");
    fprintf(file,"%u\n",points_to_faces.size());
    for (it = unique_points.begin(); it != unique_points.end(); ++it)
    {
        fprintf(file,"%lf %lf %lf\n",it->first.x,it->first.y,it->first.z);
    }
    for (uint32_t i = 0; i < points_to_faces.size(); i++)
    {
        for (uint32_t j = 0; j < points_to_faces[i].size(); j++)
        {
            fprintf(file,"%u %u\n",i,points_to_faces[i][j]);
        }
    }
    fclose(file);
}

void print_faces (std::vector<Face_Custom> faces)
{
    for (int i = 0; i < (int)faces.size(); i++)
    {
        printf("[Face %d]\n",i);
        faces[i].print();
    }
}

void usage (const char pname[])
{
  printf("%s\n",PRINT_LINE);
  printf("Usage:> %s <input_file>\n",pname);
  printf("%s\n",PRINT_LINE);
}

int main (int argc, char *argv[])
{
    if (argc-1 != 1)
    {
        usage(argv[0]);
        exit(EXIT_FAILURE);
    }

    char *mesh_filename = argv[1];
    read_faces_from_stl(mesh_filename);
    //print_faces(mesh_faces);

    write_points_from_faces_to_map();

    return 0;
}