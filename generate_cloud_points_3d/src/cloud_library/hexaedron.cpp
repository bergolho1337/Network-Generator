#include "hexaedron.h"

// Generate points inside a hexaedron with dimensions: [-side_length,side_length]x[-side_length,side_length]x[-side_length,side_length]
SET_CLOUD_GENERATOR (default_hexaedron_cloud)
{
    printf("\n[hexaedron] Generating default hexaedron cloud of points\n");

    uint32_t num_points = generator->num_points;
    std::vector<Point> *points = generator->points;
    double *rand_array = random_generator->array;
    
    double side_length;
    get_parameter_value_from_map(config->param,"side_length",&side_length);

    double pos[3];
    for (uint32_t i = 0; i < num_points; i++)
    {
        pos[0] = 2.0 * random_generator->get_value() - 1.0;
        pos[1] = 2.0 * random_generator->get_value() - 1.0;
        pos[2] = 2.0 * random_generator->get_value() - 1.0;

        // Convert to the real domain
        pos[0] *= side_length;
        pos[1] *= side_length;
        pos[2] *= side_length;

        Point p(pos[0],pos[1],pos[2]);
        points->push_back(p);
        
        printf("[hexaedron] Generating point = %u\n",i);
    }

    draw_default_hexaedron_volume(side_length);
}

void draw_default_hexaedron_volume (const double side_length)
{
    double P0[3] = {-side_length, -side_length, -side_length};
    double P1[3] = {side_length, -side_length, -side_length};
    double P2[3] = {side_length, side_length, -side_length};
    double P3[3] = {-side_length, side_length, -side_length};
    double P4[3] = {-side_length, -side_length, side_length};
    double P5[3] = {side_length, -side_length, side_length};
    double P6[3] = {side_length, side_length, side_length};
    double P7[3] = {-side_length, side_length, side_length};

    // Create the points
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    points->InsertNextPoint(P0);
    points->InsertNextPoint(P1);
    points->InsertNextPoint(P2);
    points->InsertNextPoint(P3);
    points->InsertNextPoint(P4);
    points->InsertNextPoint(P5);
    points->InsertNextPoint(P6);
    points->InsertNextPoint(P7);

    // Create a hexahedron from the points
    vtkSmartPointer<vtkHexahedron> hex = vtkSmartPointer<vtkHexahedron>::New();
    hex->GetPointIds()->SetId(0,0);
    hex->GetPointIds()->SetId(1,1);
    hex->GetPointIds()->SetId(2,2);
    hex->GetPointIds()->SetId(3,3);
    hex->GetPointIds()->SetId(4,4);
    hex->GetPointIds()->SetId(5,5);
    hex->GetPointIds()->SetId(6,6);
    hex->GetPointIds()->SetId(7,7);

    // Add the hexahedron to a cell array
    vtkSmartPointer<vtkCellArray> hexs = vtkSmartPointer<vtkCellArray>::New();
    hexs->InsertNextCell(hex);

    // Add the points and hexahedron to an unstructured grid
    vtkSmartPointer<vtkUnstructuredGrid> unstructured_grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
    unstructured_grid->SetPoints(points);
    unstructured_grid->InsertNextCell(hex->GetCellType(), hex->GetPointIds());

    vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
    writer->SetFileName("output/volume_region.vtu");
    #if VTK_MAJOR_VERSION <= 5
    writer->SetInput(unstructured_grid);
    #else
    writer->SetInputData(unstructured_grid);
    #endif
    writer->Write();
}
