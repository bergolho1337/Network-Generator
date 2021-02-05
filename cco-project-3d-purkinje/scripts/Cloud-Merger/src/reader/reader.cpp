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
        uint32_t num_cells = input_polydata->GetNumberOfCells();

        for (uint32_t i = 0; i < num_points; i++)
        {
            double pos[3];

            input_polydata->GetPoint(i,pos);

            Point point(i,pos[0],pos[1],pos[2]);

            the_points.push_back(point);
        }

        std::string array_name = "LAT";
        vtkSmartPointer<vtkFloatArray> array = vtkFloatArray::SafeDownCast(input_polydata->GetPointData()->GetArray(array_name.c_str()));
        if(array)
        {
            for(uint32_t i = 0; i < num_points; i++)
            {
                double value = array->GetValue(i);

                the_points[i].lat = value;
            }
            
        }
        array_name = "Active";
        vtkSmartPointer<vtkFloatArray> array_2 = vtkFloatArray::SafeDownCast(input_polydata->GetPointData()->GetArray(array_name.c_str()));
        if(array)
        {
            for(uint32_t i = 0; i < num_points; i++)
            {
                double value = array_2->GetValue(i);

                the_points[i].is_active = (int)value;
            }
        }
    }
    else
    {
        fprintf(stderr,"[-] ERROR! The input data must be in Polydata format!\n");
        exit(EXIT_FAILURE);
    }
}

void VTK_Reader::concatenate (VTK_Reader *input)
{
    //uint32_t offset = (this->the_points.size()+input->the_points.size()) / input->the_points.size();
    //auto it = this->the_points.begin();
    //for (uint32_t i = 0; i < input->the_points.size(); i++)
    //    this->the_points.insert(it+i*offset+offset,input->the_points[i]);
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
        double lat = this->the_points[i].lat;
        int is_active = this->the_points[i].is_active;

        printf("[%u] = (%g %g %g) {LAT=%g} {Active=%d}\n",id,x,y,z,lat,is_active);
    }
}

void VTK_Reader::write (std::string filename)
{
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    for (uint32_t i = 0; i < this->the_points.size(); i++)
    {
        double x = this->the_points[i].x;
        double y = this->the_points[i].y;
        double z = this->the_points[i].z;

        points->InsertNextPoint(x,y,z);
    }
    vtkSmartPointer<vtkCellArray> cell_array = vtkSmartPointer<vtkCellArray>::New();
    for (uint32_t i = 0; i < this->the_points.size(); i++)
    {
        vtkSmartPointer<vtkVertex> vertex = vtkSmartPointer<vtkVertex>::New();
        vertex->GetPointIds()->SetId(0,i);

        cell_array->InsertNextCell(vertex);
    }
    vtkSmartPointer<vtkFloatArray> cells_data = vtkSmartPointer<vtkFloatArray>::New();
    cells_data->SetName("LAT");
    for (uint32_t i = 0; i < this->the_points.size(); i++)
    {
        cells_data->InsertNextValue(this->the_points[i].lat);
    }
    vtkSmartPointer<vtkFloatArray> cells_data_2 = vtkSmartPointer<vtkFloatArray>::New();
    cells_data_2->SetName("Active");
    for (uint32_t i = 0; i < this->the_points.size(); i++)
    {
        cells_data_2->InsertNextValue(this->the_points[i].is_active);
    }

    vtkSmartPointer<vtkPolyData> polydata_grid = vtkSmartPointer<vtkPolyData>::New();
    polydata_grid->SetPoints(points);
    polydata_grid->SetVerts(cell_array);
    polydata_grid->GetPointData()->AddArray(cells_data);
    polydata_grid->GetPointData()->AddArray(cells_data_2);

    vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New();
    writer->SetFileName(filename.c_str());
    writer->SetInputData(polydata_grid);
    writer->Write();
}