#include "paraboloid.h"

SET_CLOUD_GENERATOR (default_paraboloid_cloud)
{
    printf("\n[sphere] Generating default paraboloid cloud of points\n");

    uint32_t num_points = generator->num_points;
    std::vector<Point> *points = generator->points;
    double *rand_array = random_generator->array;

    double a;
    get_parameter_value_from_map(config->param,"a",&a);

    double b;
    get_parameter_value_from_map(config->param,"b",&b);

    double side_length;
    get_parameter_value_from_map(config->param,"side_length",&side_length);

    double pos[3];
    for (uint32_t i = 0; i < num_points; i++)
    {
        bool point_is_inside_sphere = false;
        while (!point_is_inside_sphere)
        {
            pos[0] = 2.0 * random_generator->get_value() - 1.0;
            pos[1] = 2.0 * random_generator->get_value() - 1.0;
            pos[2] = 2.0 * random_generator->get_value() - 1.0;

            double f = a * pos[0] * pos[0] + b * pos[1] * pos[1];

            if (pos[2] > f)
                point_is_inside_sphere = true;
            else
                point_is_inside_sphere = false;
        }

        // Convert to real domain
        pos[0] *= side_length;
        pos[1] *= side_length;
        pos[2] *= side_length;
//        pos[2] -= side_length;

        Point p(pos[0],pos[1],pos[2]);
        points->push_back(p);
        
        printf("[paraboloid] Generating point = %u\n",i);
    }

    //draw_default_paraboloid_volume(a,b,side_length);
}

void draw_default_paraboloid_volume (const double a, const double b, const double side_length)
{

    // Create the quadric function definition
    vtkSmartPointer<vtkQuadric> quadric = vtkSmartPointer<vtkQuadric>::New();
    quadric->SetCoefficients(a,b,0,0,0,0,0,0,-1,0);

    vtkSmartPointer<vtkSampleFunction> sample = vtkSmartPointer<vtkSampleFunction>::New();
    sample->SetSampleDimensions(200,200,200);
    sample->SetImplicitFunction(quadric);
    double xmin = -2, xmax=2, ymin=-2, ymax=2, zmin=-1, zmax=1;
    //double xmin = -side_length, xmax=side_length, ymin=-side_length, ymax=side_length, zmin=-side_length, zmax=side_length;
    sample->SetModelBounds(xmin, xmax, ymin, ymax, zmin, zmax);

    vtkSmartPointer<vtkContourFilter> contours = vtkSmartPointer<vtkContourFilter>::New();
    contours->SetInputConnection(sample->GetOutputPort());
    contours->GenerateValues(1,1,1);

    vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->SetInputConnection(contours->GetOutputPort());
    mapper->Update();   

    // Get the reference to the Polydata
    vtkPolyData *polydata = mapper->GetInput();

    vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New();
    writer->SetFileName("output/paraboloid.vtk");
    writer->SetInputData(polydata);
    writer->Write();

    vtkSmartPointer<vtkSTLWriter> stlWriter = vtkSmartPointer<vtkSTLWriter>::New();
    stlWriter->SetFileName("output/paraboloid.stl");
    stlWriter->SetInputData(polydata);
    stlWriter->Write();
}
