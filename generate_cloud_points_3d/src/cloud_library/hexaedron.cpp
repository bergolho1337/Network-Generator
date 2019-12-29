#include "hexaedron.h"

// Generate points inside a hexaedron with dimensions: [-1,1]x[-1,1]x[-1,1]
SET_CLOUD_GENERATOR (default_hexaedron_cloud)
{
    printf("\n[hexaedron] Generating default hexaedron cloud of points\n");

    uint32_t num_points = generator->num_points;
    std::vector<Point> *points = generator->points;
    double *rand_array = random_generator->array;

    double pos[3];
    for (uint32_t i = 0; i < num_points; i++)
    {
        pos[0] = 2.0 * random_generator->get_value() - 1.0;
        pos[1] = 2.0 * random_generator->get_value() - 1.0;
        pos[2] = 2.0 * random_generator->get_value() - 1.0;

        Point p(pos[0],pos[1],pos[2]);
        points->push_back(p);
        
        printf("[hexaedron] Generating point = %u\n",i);
    }

    draw_default_hexaedron_volume(1.0);
}

void draw_default_hexaedron_volume (const double side_length)
{
    double P0[3] = {-1.0, -1.0, -1.0};
    double P1[3] = {1.0, -1.0, -1.0};
    double P2[3] = {1.0, 1.0, -1.0};
    double P3[3] = {-1.0, 1.0, -1.0};
    double P4[3] = {-1.0, -1.0, 1.0};
    double P5[3] = {1.0, -1.0, 1.0};
    double P6[3] = {1.0, 1.0, 1.0};
    double P7[3] = {-1.0, 1.0, 1.0};

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
    writer->SetFileName("output/hexaedron_volume_region.vtu");
    #if VTK_MAJOR_VERSION <= 5
    writer->SetInput(unstructured_grid);
    #else
    writer->SetInputData(unstructured_grid);
    #endif
    writer->Write();
}
