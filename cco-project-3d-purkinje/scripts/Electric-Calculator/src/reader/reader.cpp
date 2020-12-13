#include "reader.h"

VTK_Reader::VTK_Reader (std::string filename)
{
    // Read all the data from the VTK file
    vtkSmartPointer<vtkGenericDataObjectReader> reader = vtkSmartPointer<vtkGenericDataObjectReader>::New();
    reader->SetFileName(filename.c_str());
    reader->Update();

    if(reader->IsFilePolyData())
    {
        vtkSmartPointer<vtkPolyData> input_polydata = reader->GetPolyDataOutput();

        uint32_t num_points = input_polydata->GetNumberOfPoints();
        uint32_t num_cells = input_polydata->GetNumberOfCells();

        for (uint32_t i = 0; i < num_points; i++)
        {
            double pos[3];

            input_polydata->GetPoint(i,pos);

            Point p(i,pos[0],pos[1],pos[2]);

            this->the_points.push_back(p);

            //printf("[%u] = (%g %g %g)\n",i,pos[0],pos[1],pos[2]);
        }

        for (uint32_t i = 0; i < num_cells; i++)
        {
            vtkCell *cell = input_polydata->GetCell(i);

            vtkLine *line = dynamic_cast<vtkLine*>(cell);

            uint32_t src = line->GetPointIds()->GetId(0);
            uint32_t dest = line->GetPointIds()->GetId(1);

            Line l(src,dest);

            this->the_lines.push_back(l);

            //printf("[%u] %u --> %u\n",i,src,dest);
        }
    }
}

vtkSmartPointer<vtkPolyData> VTK_Reader::convert_to_polydata ()
{
    vtkSmartPointer<vtkPolyData> output_polydata = vtkSmartPointer<vtkPolyData>::New();   
    
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    for (uint32_t i = 0; i < this->the_points.size(); i++)
    {
        double x = this->the_points[i].x;
        double y = this->the_points[i].y;
        double z = this->the_points[i].z;

        points->InsertNextPoint(x,y,z);
    }

    vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();
    for (uint32_t i = 0; i < this->the_lines.size(); i++)
    {
        uint32_t src = this->the_lines[i].src;
        uint32_t dest = this->the_lines[i].dest;

        vtkSmartPointer<vtkLine> line  = vtkSmartPointer<vtkLine>::New();
        line->GetPointIds()->SetId(0,src);
        line->GetPointIds()->SetId(1,dest);

        lines->InsertNextCell(line);
    }

    output_polydata->SetPoints(points);
    output_polydata->SetLines(lines);

    return output_polydata;
} 

uint32_t VTK_Reader::search_position (const double pos[])
{
    uint32_t closest_index = 0;
    double closest_dist = __DBL_MAX__;
    for (uint32_t i = 0; i < this->the_points.size(); i++)
    {
        double x = this->the_points[i].x;
        double y = this->the_points[i].y;
        double z = this->the_points[i].z;

        double dist = sqrt( pow(x-pos[0],2) + pow(y-pos[1],2) + pow(z-pos[2],2) );
        if (dist < closest_dist)
        {
            closest_dist = dist;
            closest_index = i;
        }
    }
    return closest_index;
}

void read_points (std::string filename, std::vector<Point> &the_points)
{
    uint32_t np;
    FILE *file = fopen(filename.c_str(),"r");
    fscanf(file,"%u",&np);
    for (uint32_t i = 0 ; i < np; i++)
    {
        double x, y, z;
        fscanf(file,"%lf %lf %lf",&x,&y,&z);
        Point p(i,x,y,z);
        the_points.push_back(p);
    }
    fclose(file);
}

void VTK_Reader::print ()
{
    printf("======= POINTS =======\n");
    for (uint32_t i = 0; i < this->the_points.size(); i++)
    {
        uint32_t id = this->the_points[i].id;
        double x = this->the_points[i].x;
        double y = this->the_points[i].y;
        double z = this->the_points[i].z;

        printf("[%u] = (%g %g %g)\n",id,x,y,z);
    }

    printf("======= LINES =======\n");
    for (uint32_t i = 0; i < this->the_lines.size(); i++)
    {
        uint32_t id = i;
        uint32_t src = this->the_lines[i].src;
        uint32_t dest = this->the_lines[i].dest;

        printf("[%u] %u --> %u\n",i,src,dest);
    }
}

void VTK_Reader::write (std::string filename)
{
    vtkSmartPointer<vtkPolyData> output_polydata = convert_to_polydata();

    vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New();
    writer->SetFileName(filename.c_str());
    writer->SetInputData(output_polydata);
    writer->Write();
}
