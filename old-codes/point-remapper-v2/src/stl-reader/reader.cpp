#include "reader.h"

uint32_t counter_points_table = 0;

void insert_point_into_table (map<Point_Table,uint32_t> &point_table, const double pos[])
{
    map<Point_Table,uint32_t>::iterator it;
    Point_Table point(pos[0],pos[1],pos[2]);

    it = point_table.find(point);
    if (it == point_table.end())
    {
        point_table.insert(pair<Point_Table,uint32_t>(point,counter_points_table));
        counter_points_table++;
    }
}

void read_face (FILE *file, vector<Face> &faces, map<Point_Table,uint32_t> &point_table)
{
    map<Point_Table,uint32_t>::iterator it;
    char str[200];
    double n[3];
    double a[3], b[3], c[3];
    Point *v1, *v2, *v3;

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

    insert_point_into_table(point_table,a);
    insert_point_into_table(point_table,b);
    insert_point_into_table(point_table,c);

    Point_Table p1(a[0],a[1],a[2]);
    it = point_table.find(p1);
    if (it != point_table.end())    // Point was found in the table
    {
        uint32_t point_id = it->second;
        v1 = new Point(point_id,a[0],a[1],a[2]);
    }
    Point_Table p2(b[0],b[1],b[2]);
    it = point_table.find(p2);
    if (it != point_table.end())    // Point was found in the table
    {
        uint32_t point_id = it->second;
        v2 = new Point(point_id,b[0],b[1],b[2]);
    }
    Point_Table p3(c[0],c[1],c[2]);
    it = point_table.find(p3);
    if (it != point_table.end())    // Point was found in the table
    {
        uint32_t point_id = it->second;
        v3 = new Point(point_id,c[0],c[1],c[2]);
    }

    Face new_face(v1,v2,v3,n);

    faces.push_back(new_face);
}

void read_faces (const char *filename, vector<Face> &faces, map<Point_Table,uint32_t> &point_table)
{
    char str[200];
    FILE *file = fopen(filename,"r");

    while (fscanf(file,"%s",str) != EOF)
    {
        if (strcmp(str,"facet") == 0)
            read_face(file,faces,point_table);
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