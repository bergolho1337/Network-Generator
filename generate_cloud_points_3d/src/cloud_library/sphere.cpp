#include "sphere.h"

// Generate points inside a sphere of radius 1, which is centered on (0,0,0)
SET_CLOUD_GENERATOR (default_sphere_cloud)
{
    printf("\n[sphere] Generating default sphere cloud of points\n");

    uint32_t num_points = generator->num_points;
    std::vector<Point> *points = generator->points;
    double *rand_array = random_generator->array;

    double radius;
    get_parameter_value_from_map(config->param,"radius",&radius);

    double pos[3];
    for (uint32_t i = 0; i < num_points; i++)
    {
        bool point_is_inside_sphere = false;
        while (!point_is_inside_sphere)
        {
            pos[0] = 2.0 * random_generator->get_value() - 1.0;
            pos[1] = 2.0 * random_generator->get_value() - 1.0;
            pos[2] = 2.0 * random_generator->get_value() - 1.0;

            // Convert to the real domain
            pos[0] *= radius;
            pos[1] *= radius;
            pos[2] *= radius;
            
            if (sqrt(pos[0] * pos[0] + pos[1] * pos[1] + pos[2] * pos[2]) > radius)
                point_is_inside_sphere = false;
            else
                point_is_inside_sphere = true;
        }

        Point p(pos[0],pos[1],pos[2]);
        points->push_back(p);
        
        printf("[sphere] Generating point = %u\n",i);
    }

    draw_default_sphere_volume(radius);
}

void draw_default_sphere_volume (const double radius)
{
    vtkSmartPointer<vtkSphereSource> sphereSource = vtkSmartPointer<vtkSphereSource>::New();
    
    sphereSource->SetRadius(radius);
    sphereSource->SetCenter(0,0,0);
    sphereSource->SetPhiResolution(100);
    sphereSource->SetThetaResolution(100);

    vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
    writer->SetFileName("output/volume_region.vtp");
    writer->SetInputConnection(sphereSource->GetOutputPort());
    writer->Write();
}
