#include "custom.h"

SET_WALKER_MOVE_FUNCTION (move)
{
    uint32_t src_face_index, dest_face_index;
    std::map<Point_Custom,uint32_t>::iterator it, it2;
    Point_Custom tmp(the_walker->pos[0],the_walker->pos[1],the_walker->pos[2]);

    // Find the point which is at the Walker current position in the map
    it = unique_points.find(tmp);
    if (it == unique_points.end())
    {
        fprintf(stderr,"[-] ERROR! Point not found!\n");
        exit(EXIT_FAILURE);
    }

    // Get the face which contains point and sort a neighbour face
    src_face_index = it->second;
    uint32_t num_neighbour_faces = points_to_faces[src_face_index].size();
    dest_face_index = points_to_faces[src_face_index][rand() % num_neighbour_faces];

    // Next, sort a vertex from the destination face to move to
    uint32_t vertex_index = rand() % 3;

    switch (vertex_index)
    {
        // Vertex 1
        case 0: {
                    the_walker->pos[0] = mesh_faces[dest_face_index].v1->x;
                    the_walker->pos[1] = mesh_faces[dest_face_index].v1->y;
                    the_walker->pos[2] = mesh_faces[dest_face_index].v1->z;
                    break;
                }
        // Vertex 2
        case 1: {
                    the_walker->pos[0] = mesh_faces[dest_face_index].v2->x;
                    the_walker->pos[1] = mesh_faces[dest_face_index].v2->y;
                    the_walker->pos[2] = mesh_faces[dest_face_index].v2->z;
                    break;
                }
        // Vertex 3
        case 2: {
                    the_walker->pos[0] = mesh_faces[dest_face_index].v3->x;
                    the_walker->pos[1] = mesh_faces[dest_face_index].v3->y;
                    the_walker->pos[2] = mesh_faces[dest_face_index].v3->z;
                    break;
                }
    }
}

SET_WALKER_RESPAWN_FUNCTION (respawn)
{
    if (first_call)
    {
        char *mesh_filename = the_walker_config->mesh_filename;
        char *map_filename = the_walker_config->map_filename;

        printf("[custom] Reading STL file with the mesh nodes ...\n");
        read_faces_from_stl(mesh_filename);

        if (!map_filename)
        {
            printf("[custom] Calculating the points to faces mapping ...\n");
            insert_points_from_faces_to_map();
        }
        else
        {
            printf("[custom] Reading the points to faces mapping from file '%s' ...\n",map_filename);
            read_points_from_faces_to_map(map_filename);
        }

        first_call = false;
    }

    uint32_t face_index = rand() % mesh_faces.size();
    uint32_t vertex_index = rand() % 3;
    switch (vertex_index)
    {
        // Vertex 1
        case 0: {
                    pos[0] = mesh_faces[face_index].v1->x;
                    pos[1] = mesh_faces[face_index].v1->y;
                    pos[2] = mesh_faces[face_index].v1->z;
                    break;
                }
        // Vertex 2
        case 1: {
                    pos[0] = mesh_faces[face_index].v2->x;
                    pos[1] = mesh_faces[face_index].v2->y;
                    pos[2] = mesh_faces[face_index].v2->z;
                    break;
                }
        // Vertex 3
        case 2: {
                    pos[0] = mesh_faces[face_index].v3->x;
                    pos[1] = mesh_faces[face_index].v3->y;
                    pos[2] = mesh_faces[face_index].v3->z;
                    break;
                }
    }
}

SET_WALKER_DOMAIN_DRAW_FUNCTION (draw)
{
    
}

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

void print_faces (std::vector<Face_Custom> faces)
{
    for (int i = 0; i < (int)faces.size(); i++)
    {
        printf("[Face %d]\n",i);
        faces[i].print();
    }
}

void insert_points_from_faces_to_map ()
{
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
    
    // DEBUG
/*
    for (uint32_t i = 0; i < points_to_faces.size(); i++)
    {
        printf("Point %u -- ",i);
        for (uint32_t j = 0; j < points_to_faces[i].size(); j++)
        {
            printf("%u ",points_to_faces[i][j]);
        }
        printf("\n");
    }
*/
}

void read_points_from_faces_to_map (const char filename[])
{
    FILE *file = fopen(filename,"r");
    if (!file)
    {
        fprintf(stderr,"[custom] ERROR! Could not open file '%s'!\n",filename);
        exit(EXIT_FAILURE);
    }

    // Read the number of points
    uint32_t num_points;
    fscanf(file,"%u",&num_points);

    // Insert the unique points
    uint32_t unique_num_points = 0;
    for (uint32_t i = 0; i < mesh_faces.size(); i++)
    {
        Point_Custom *v1 = mesh_faces[i].v1;
        Point_Custom *v2 = mesh_faces[i].v2;
        Point_Custom *v3 = mesh_faces[i].v3;

        auto it = unique_points.find(*v1);
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

    // Initialize the map 
    points_to_faces.assign(unique_points.size(),std::vector<uint32_t>());

    // Read all the links
    uint32_t point_index, face_index;
    while (fscanf(file,"%u %u",&point_index,&face_index) != EOF)
    {
        points_to_faces[point_index].push_back(face_index);
    }

    fclose(file);
}