#include "tree.h"

Tree::Tree ()
{

}

void Tree::grow ()
{
    // Make the root
    Walker *root = new Walker(WIDTH/2.0,HEIGHT/2.0,0.0);
    this->the_tree.push_back(root);

    // Add the Walkers
    std::vector<Walker*> the_others;
    for (uint32_t i = 0; i < MAX_NUMBER_OF_WALKERS; i++)
    {
        Walker *walker = new Walker();
        //this->the_tree.push_back(walker);
        the_others.push_back(walker);
    }

    // Main iteration loop 
    for (uint32_t iter = 0; iter < MAX_NUMBER_OF_ITERATIONS; iter++)
    {
        print_progress_bar(iter,MAX_NUMBER_OF_ITERATIONS);

        // Move each Walker using the Random Walk
        for (uint32_t i = the_others.size()-1; i > 0; i--)
        {
            //printf("%u\n",i);
            the_others[i]->walk();

            uint32_t stuck_index = the_others[i]->is_stuck(this->the_tree);
            if (stuck_index != this->the_tree.size())
            {
                uint32_t new_index = this->the_tree.size();
                Segment *segment = new Segment(stuck_index,new_index);
                this->the_segments.push_back(segment);

                this->the_tree.push_back(the_others[i]);
                the_others.erase(the_others.begin()+i);

                write_to_vtk(this->the_tree.size());
            }
        }
        printf("\n");

        // Add more walkers as the tree grows ...
        while (the_others.size() < MAX_NUMBER_OF_WALKERS)
        {
            Walker *walker = new Walker();
            the_others.push_back(walker);
        }
    }

}

void Tree::write_to_vtk (const uint32_t iter)
{
    uint32_t num_points = this->the_tree.size();
    uint32_t num_lines = this->the_segments.size();

    char filename[50];
    sprintf(filename,"output/tree_%u.vtk",iter);
    FILE *file = fopen(filename,"w+");

    fprintf(file,"# vtk DataFile Version 3.0\n");
    fprintf(file,"Tree\n");
    fprintf(file,"ASCII\n");
    fprintf(file,"DATASET POLYDATA\n");
    fprintf(file,"POINTS %u float\n",num_points);
    
    // Write the points
    for (uint32_t i = 0; i < this->the_tree.size(); i++)    
    {
        fprintf(file,"%g %g %g\n",this->the_tree[i]->pos[0],this->the_tree[i]->pos[1],this->the_tree[i]->pos[2]);
    }

    fprintf(file,"LINES %u %u\n",num_lines,num_lines*3);

    // Write the lines
    for (uint32_t i = 0; i < num_lines; i++)
    {
        fprintf(file,"2 %u %u\n",this->the_segments[i]->src,this->the_segments[i]->dest);
    }
    
    fclose(file);
}

void Tree::write_to_vtp ()
{
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();       
    vtkSmartPointer<vtkCellArray> cell_array = vtkSmartPointer<vtkCellArray>::New();
        
    uint32_t counter = 0;
    for (uint32_t i = 0; i < the_tree.size(); i++)
    {
        double *pos = the_tree[i]->pos;

        points->InsertNextPoint(pos[0]-12,pos[1]-12,pos[2]);
        points->InsertNextPoint(pos[0]+12,pos[1]-12,pos[2]);
        points->InsertNextPoint(pos[0]-12,pos[1]-12,pos[2]+12);
        points->InsertNextPoint(pos[0]+12,pos[1]-12,pos[2]+12);
        points->InsertNextPoint(pos[0]+12,pos[1]+12,pos[2]);
        points->InsertNextPoint(pos[0]-12,pos[1]+12,pos[2]);
        points->InsertNextPoint(pos[0]+12,pos[1]+12,pos[2]+12);
        points->InsertNextPoint(pos[0]-12,pos[1]+12,pos[2]+12);

        vtkSmartPointer<vtkHexahedron> hexahedron_1 = vtkSmartPointer<vtkHexahedron>::New();
        hexahedron_1->GetPointIds()->SetId(0,counter+0);
        hexahedron_1->GetPointIds()->SetId(1,counter+1);
        hexahedron_1->GetPointIds()->SetId(2,counter+2);
        hexahedron_1->GetPointIds()->SetId(3,counter+3);
        hexahedron_1->GetPointIds()->SetId(4,counter+4);
        hexahedron_1->GetPointIds()->SetId(5,counter+6);
        hexahedron_1->GetPointIds()->SetId(6,counter+7);
        hexahedron_1->GetPointIds()->SetId(7,counter+8);

        cell_array->InsertNextCell(hexahedron_1);

        counter += 8;
   
    }

    vtkSmartPointer<vtkUnstructuredGrid> unstructured_grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
    unstructured_grid->SetPoints(points);
    unstructured_grid->SetCells(VTK_HEXAHEDRON,cell_array);

    // Write the polydata to a file
    vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
    writer->SetFileName("output/tree.vtu");
    writer->SetInputData(unstructured_grid);
    writer->Write();

}