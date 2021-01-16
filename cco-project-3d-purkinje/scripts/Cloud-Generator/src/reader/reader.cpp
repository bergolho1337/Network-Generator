#include "reader.h"

VTK_Reader::VTK_Reader (std::string filename)
{
    vtkSmartPointer<vtkGenericDataObjectReader> reader = vtkSmartPointer<vtkGenericDataObjectReader>::New();
    reader->SetFileName(filename.c_str());
    reader->Update();

    if(reader->IsFilePolyData())
    {
        vtkSmartPointer<vtkPolyData> input_polydata = reader->GetPolyDataOutput();

        uint32_t num_points = input_polydata->GetNumberOfPoints();
       
        for (uint32_t i = 0; i < num_points; i++)
        {
            double pos[3];

            input_polydata->GetPoint(i,pos);

            Point p(i,pos[0],pos[1],pos[2]);

            this->the_points.push_back(p);
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

    vtkSmartPointer<vtkCellArray> vertexes = vtkSmartPointer<vtkCellArray>::New();
    for (uint32_t i = 0; i < this->the_points.size(); i++)
    {
        vtkSmartPointer<vtkVertex> vertex  = vtkSmartPointer<vtkVertex>::New();
        vertex->GetPointIds()->SetId(0,i);

        vertexes->InsertNextCell(vertex);
    }

    output_polydata->SetPoints(points);
    output_polydata->SetVerts(vertexes);

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

void VTK_Reader::concatenate (VTK_Reader *input)
{
    for (uint32_t i = 0; i < input->the_points.size(); i++)
        this->the_points.push_back(input->the_points[i]);
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
}

void VTK_Reader::write (std::string filename)
{
    vtkSmartPointer<vtkPolyData> output_polydata = convert_to_polydata();

    vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New();
    writer->SetFileName(filename.c_str());
    writer->SetInputData(output_polydata);
    writer->Write();
}

bool VTK_Reader::check_duplicates ()
{
    for (uint32_t i = 0; i < this->the_points.size(); i++)
    {
        double p1[3];
        p1[0] = this->the_points[i].x;
        p1[1] = this->the_points[i].y;
        p1[2] = this->the_points[i].z;
        for (uint32_t j = i+1; j < this->the_points.size(); j++)
        {
            double p2[3];
            p2[0] = this->the_points[j].x;
            p2[1] = this->the_points[j].y;
            p2[2] = this->the_points[j].z;
            double dist = sqrt(powf(p1[0]-p2[0],2)+powf(p1[1]-p2[1],2)+powf(p1[2]-p2[2],2));
            if (dist < 1.0e-08)
            {
                printf("Duplicate points: %u <--> %u\n",i,j);
                return true;
            }
        }
    }
    return false;
}