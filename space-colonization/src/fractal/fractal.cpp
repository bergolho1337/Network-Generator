#include "fractal.h"

Fractal_Tree::Fractal_Tree ()
{
    this->max_iterations = -1;

    this->initial_angle = 0.0;
    this->initial_length = -1.0;
    this->length_decrease_ratio = -1.0;
    this->angle_decrease_ratio = -1.0;
    
    for (uint32_t i = 0; i < 3; i++)
        this->root_pos[i] = 0.0;
}

void Fractal_Tree::make_root (std::queue<Line> &q)
{
    printf("[fractal] Making root at position --> (%g,%g,%g)\n",this->root_pos[0],this->root_pos[1],this->root_pos[2]);

    double *u = this->root_pos;
    double initial_length = this->initial_length;
    double initial_diameter = this->initial_diameter;

    // Direction vector
    double d[3] = {1,0,0};

    // Root source point
    Point root_src(0,u[0],u[1],u[2]);

    // Root destination point
    double v[3];
    for (uint32_t i = 0; i < 3; i++)
        v[i] = u[i] + d[i]*initial_length;
    Point root_dest(1,v[0],v[1],v[2]);

    // Insert the points into the Point array
    this->the_points.push_back(root_src);
    this->the_points.push_back(root_dest);

    // Create and insert the root segment into the Line array
    Line root_line(0,0,1,initial_diameter);
    this->the_lines.push_back(root_line);

    // Enqueue the root segment
    q.push(root_line);
}

void Fractal_Tree::grow_network ()
{
    uint32_t max_iterations = this->max_iterations;
    double initial_length = this->initial_length;
    double initial_angle = this->initial_angle;
    double initial_diameter = this->initial_diameter;
    double length_ratio = this->length_decrease_ratio;
    double angle_ratio = this->angle_decrease_ratio;
    double diameter_ratio = this->diameter_decrease_ratio;

    // Growing queue of segments
    std::queue<Line> q;

    // Build the root segment
    make_root(q);

    // MAIN LOOP
    for (uint32_t iter = 1; iter <= max_iterations; iter++)
    {
        uint32_t num_segments_to_grow = q.size();

        double iteration_length = initial_length * powf(length_ratio,iter);
        double iteration_angle = initial_angle * powf(angle_ratio,iter);
        double iteration_diameter = initial_diameter * powf(diameter_ratio,iter);
        printf("[fractal] Iteration %u -- Length = %g -- Angle = %g -- Diameter = %g\n",iter,iteration_length,iteration_angle,iteration_diameter);

        for (uint32_t k = 0; k < num_segments_to_grow; k++)
        {
            Line cur_seg = q.front(); q.pop();
            uint32_t src_id = cur_seg.src;
            uint32_t dest_id = cur_seg.dest;

            // Copy the positions from the Points that define the current Line
            double u[3], v[3];
            u[0] = this->the_points[src_id].x;
            u[1] = this->the_points[src_id].y;
            u[2] = this->the_points[src_id].z;

            v[0] = this->the_points[dest_id].x;
            v[1] = this->the_points[dest_id].y;
            v[2] = this->the_points[dest_id].z;

            // Calculate the normal growth direction
            //double d[3];
            //calculate_unitary_vector(u,v,d);
            double d[3] = {1,0,0};  // Constant direction vector

            // Calculate the positions from the two new offsprings
            double b1[3], b2[3];
            calculate_new_branch_position(b1,v,d,iteration_length,iteration_angle);
            calculate_new_branch_position(b2,v,d,iteration_length,-iteration_angle);

            // Build the new points
            uint32_t p1_id = this->the_points.size();
            uint32_t p2_id = this->the_points.size()+1;
            Point p1(p1_id,b1[0],b1[1],b1[2]);
            Point p2(p2_id,b2[0],b2[1],b2[2]);

            // Insert then into the Points array
            this->the_points.push_back(p1);
            this->the_points.push_back(p2);

            // Build the new branches
            uint32_t l1_id = this->the_lines.size();
            uint32_t l2_id = this->the_lines.size()+1;
            Line l1(l1_id,dest_id,p1_id,iteration_diameter);
            Line l2(l2_id,dest_id,p2_id,iteration_diameter);

            // Insert then into the Lines array
            this->the_lines.push_back(l1);
            this->the_lines.push_back(l2);

            // Enqueue the new branches for the next growth iteration
            q.push(l1);
            q.push(l2);
        }
    }
}

void Fractal_Tree::print ()
{
    printf("=============== POINTS ================\n");
    for (uint32_t i = 0; i < this->the_points.size(); i++)
        this->the_points[i].print();
    printf("=============== LINES ================\n");
    for (uint32_t i = 0; i < this->the_lines.size(); i++)
        this->the_lines[i].print();
}

void Fractal_Tree::write ()
{
    FILE *file = fopen("outputs/fractal_tree.vtk","w+");

    fprintf(file,"# vtk DataFile Version 3.0\n");
    fprintf(file,"Tree\n");
    fprintf(file,"ASCII\n");
    fprintf(file,"DATASET POLYDATA\n");

    fprintf(file,"POINTS %u float\n",this->the_points.size());
    for (uint32_t i = 0; i < this->the_points.size(); i++)
        fprintf(file,"%g %g %g\n",this->the_points[i].x,this->the_points[i].y,this->the_points[i].z);

    fprintf(file,"LINES %u %u\n",this->the_lines.size(),this->the_lines.size()*3);
    for (uint32_t i = 0; i < this->the_lines.size(); i++)
        fprintf(file,"2 %u %u\n",this->the_lines[i].src,this->the_lines[i].dest);
    
    fprintf(file,"CELL_DATA %u\n",this->the_lines.size());
    fprintf(file,"SCALARS diameter float\n");
    fprintf(file,"LOOKUP_TABLE default\n");
    for (uint32_t i = 0; i < this->the_lines.size(); i++)
        fprintf(file,"%g\n",this->the_lines[i].diameter);
    
    fclose(file);
}