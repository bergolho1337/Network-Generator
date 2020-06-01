#include "tree.h"

Tree::Tree (User_Options *options)
{
    this->max_number_iterations = options->number_of_iterations;
    memcpy(root_pos,options->root_pos,sizeof(double)*3);
    memcpy(root_dir,options->root_dir_pos,sizeof(double)*3);

    char *cloud_points_filename = options->cloud_points_filename;
    if (cloud_points_filename)
    {
        read_cloud_points(cloud_points_filename);
    }
    else
    {
        fprintf(stderr,"[ERROR] Reading cloud of points filename '%s'\n",cloud_points_filename);
        exit(EXIT_FAILURE);
    }

    uint32_t leaves_offset = options->leaves_offset;
    create_leaves(leaves_offset);

    make_root();
}

void Tree::create_leaves (const double leaves_offset)
{
    uint32_t num_points = this->the_cloud_points.size();
    for (uint32_t i = 0; i < num_points; i = i + leaves_offset)
    {
        Leaf leaf(this->the_cloud_points[i].pos[0],this->the_cloud_points[i].pos[1],this->the_cloud_points[i].pos[2]);

        this->the_leaves.push_back(leaf);
    }

    // DEBUG
    //write_leaves(0);
}

void Tree::make_root ()
{
    // Create the root segment
    double *root_pos = this->root_pos;
    double *root_dir = this->root_dir;
    Branch root(0,root_pos[0],root_pos[1],root_pos[2],\
                root_dir[0],root_dir[1],root_dir[2],\
                -1);
    this->the_branches.push_back(root);

    Branch *current = &this->the_branches[0];
    uint32_t counter = 0;
    bool found = false;

    // Until we have not found a close enough leaf to the root branch, so we could start the algorithm, 
    // we keep growing the root segment on its initial direction
    while (!found)
    {
        for (uint32_t i = 0; i < this->the_leaves.size(); i++)
        {
            double dist = calculate_distance(current->pos[0],current->pos[1],current->pos[2],\
                                            this->the_leaves[i].pos[0],this->the_leaves[i].pos[1],this->the_leaves[i].pos[2]);
            
            if (dist < MAX_DISTANCE)
            {
                printf("Found --> (%g,%g,%g) || distance = %g\n",this->the_leaves[i].pos[0],this->the_leaves[i].pos[1],this->the_leaves[i].pos[2],dist);
                found = true;
            }
        }

        // If we not found a feasible leaf them keep growing the root segment
        if (!found)
        {
            double new_pos[3];
            current->get_next_branch_position(new_pos);
            Branch branch(counter+1,new_pos[0],new_pos[1],new_pos[2],\
                            current->dir[0],current->dir[1],current->dir[2],\
                            current->id);
            this->the_branches.push_back(branch);

            current = &branch;
            counter++;
        }
    }

    write_branches(0);
}

void Tree::grow_network ()
{
    // Figure out which is the closest branch from each leaf
    for (uint32_t i = 0; i < this->the_leaves.size(); i++)
    {
        Leaf leaf = this->the_leaves[i];

        int closest = -1;
        double closest_dist = -1;
        for (uint32_t j = 0; j < this->the_branches.size(); j++)
        {
            Branch branch = this->the_branches[j];
            
            double dist = calculate_distance(leaf.pos[0],leaf.pos[1],leaf.pos[2],branch.pos[0],branch.pos[1],branch.pos[2]);

            // The leaf has been reached and it is done. Mark it for removal
            if (dist < MIN_DISTANCE)
            {
                this->the_leaves[i].is_reached = true;
                closest = -1;
                break;
            }
            else if (dist > MAX_DISTANCE)
            {   
                // Do nothing
            }
            // Update the closest branch since (dist >= MIN_DISTANCE && dist <= MAX_DISTANCE)
            else if (closest == -1 || dist < closest_dist)
            {
                closest = branch.id;
                closest_dist = dist;
                
                //printf("Leaf %u -- Closest branch = %u -- Distance = %g\n",i,closest->id,min_dist);
            }
        }

        // Create a new branch based on the closest branch of the current leaf
        if (closest != -1)
        {
            //printf("Leaf %u -- Closest branch = %u\n",i,closest->id);
            double new_dir[3];
            calculate_difference_vector(leaf.pos,this->the_branches[closest].pos,new_dir);
            normalize_vector(new_dir);

            for (uint32_t k = 0; k < 3; k++)
                this->the_branches[closest].dir[k] += new_dir[k];
            this->the_branches[closest].counter++;
        }
    }

    // Delete any leaf that was reached
    uint32_t total_leaves = this->the_leaves.size();
    for (int i = total_leaves-1; i >= 0; i--)
    {
        Leaf leaf = this->the_leaves[i];
        if (leaf.is_reached)
            this->the_leaves.erase(this->the_leaves.begin()+i);
    }
    
    // Creating new branches
    uint32_t total_branches = this->the_branches.size();
    for (int i = total_branches-1; i >= 0; i--)
    {
        if (this->the_branches[i].counter > 0)
        {
            //printf("Branch %u -- counter = %u\n",i,this->the_branches[i].counter);
            
            // Average the direction vector by the number of leaves that particular branch is close by
            for (uint32_t j = 0; j < 3; j++)
                this->the_branches[i].dir[j] /= (this->the_branches[i].counter+1);
            
            // Create a new branch from the previous one
            double new_pos[3];
            this->the_branches[i].get_next_branch_position(new_pos);

            Branch new_branch(this->the_branches.size(),new_pos[0],new_pos[1],new_pos[2],\
                            this->the_branches[i].dir[0],this->the_branches[i].dir[1],this->the_branches[i].dir[2],\
                            this->the_branches[i].id);
            
            // Add this branch to the vector
            this->the_branches.push_back(new_branch);

            // Remember to reset the branch direction and counter for the next iteration
            this->the_branches[i].reset();

        }
    }    
    //printf("\n");
}

void Tree::read_cloud_points (const char filename[])
{
    double pos[3];
    uint32_t num_points;

    FILE *file = fopen(filename,"r");
    fscanf(file,"%u",&num_points);

    for (uint32_t i = 0; i < num_points; i++)
    {
        fscanf(file,"%lf %lf %lf",&pos[0],&pos[1],&pos[2]);

        Point p(i,pos[0],pos[1],pos[2]);

        this->the_cloud_points.push_back(p);
    }   

    fclose(file);

}

void Tree::generate ()
{
    uint32_t max_number_iterations = this->max_number_iterations;

    // MAIN LOOP
    for (uint32_t k = 0; k < max_number_iterations; k++)
    {
        printf("[generate] Growing iteration %u ...\n",k);
		grow_network();

        // DEBUG
		//write_leaves(k);
        //write_branches(k);
    }

    write_branches(0);
}

void Tree::write_branches (uint32_t iter)
{
    char filename[200];
    sprintf(filename,"outputs/space_colonization_tree_%u.vtk",iter);
    FILE *file = fopen(filename,"w+");

    fprintf(file,"# vtk DataFile Version 3.0\n");
    fprintf(file,"Tree\n");
    fprintf(file,"ASCII\n");
    fprintf(file,"DATASET POLYDATA\n");

    fprintf(file,"POINTS %u float\n",this->the_branches.size());
    for (uint32_t i = 0; i < this->the_branches.size(); i++)
        fprintf(file,"%g %g %g\n",this->the_branches[i].pos[0],this->the_branches[i].pos[1],this->the_branches[i].pos[2]);

    fprintf(file,"LINES %u %u\n",this->the_branches.size()-1,(this->the_branches.size()-1)*3);
    for (uint32_t i = 0; i < this->the_branches.size(); i++)
    {
        if (this->the_branches[i].parent != -1)
        {
            uint32_t src_id = this->the_branches[i].parent;
            uint32_t dest_id = this->the_branches[i].id;
            fprintf(file,"2 %u %u\n",src_id,dest_id);
        }
    }
    
    fclose(file);
}

void Tree::write_leaves (uint32_t iter)
{
    char filename[200];
    sprintf(filename,"outputs/leaves_%u.vtk",iter);
    FILE *file = fopen(filename,"w+");

    fprintf(file,"# vtk DataFile Version 3.0\n");
    fprintf(file,"Tree\n");
    fprintf(file,"ASCII\n");
    fprintf(file,"DATASET POLYDATA\n");

    fprintf(file,"POINTS %u float\n",this->the_leaves.size());
    for (uint32_t i = 0; i < this->the_leaves.size(); i++)
        fprintf(file,"%g %g %g\n",this->the_leaves[i].pos[0],this->the_leaves[i].pos[1],this->the_leaves[i].pos[2]);

    fprintf(file,"VERTICES %u %u\n",this->the_leaves.size(),this->the_leaves.size()*2);
    for (uint32_t i = 0; i < this->the_leaves.size(); i++)
    {
        fprintf(file,"1 %u\n",i);
    }
    
    fclose(file);
}

void draw_canvas ()
{
    uint32_t width = WIDTH;
    uint32_t height = HEIGHT;

    vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
    
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    points->InsertNextPoint(0,height,0);
    points->InsertNextPoint(width,height,0);
    points->InsertNextPoint(width,0,0);
    points->InsertNextPoint(0,0,0);
    polydata->SetPoints(points);

    vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();
    vtkSmartPointer<vtkLine> line1 = vtkSmartPointer<vtkLine>::New();
    line1->GetPointIds()->SetId(0,0);
    line1->GetPointIds()->SetId(1,1);
    vtkSmartPointer<vtkLine> line2 = vtkSmartPointer<vtkLine>::New();
    line2->GetPointIds()->SetId(0,1);
    line2->GetPointIds()->SetId(1,2);
    vtkSmartPointer<vtkLine> line3 = vtkSmartPointer<vtkLine>::New();
    line3->GetPointIds()->SetId(0,2);
    line3->GetPointIds()->SetId(1,3);
    vtkSmartPointer<vtkLine> line4 = vtkSmartPointer<vtkLine>::New();
    line4->GetPointIds()->SetId(0,3);
    line4->GetPointIds()->SetId(1,0);
    lines->InsertNextCell(line1);
    lines->InsertNextCell(line2);
    lines->InsertNextCell(line3);
    lines->InsertNextCell(line4);
    polydata->SetLines(lines);

    vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
    writer->SetFileName("outputs/canvas.vtp");
    writer->SetInputData(polydata);
    writer->Write();
}