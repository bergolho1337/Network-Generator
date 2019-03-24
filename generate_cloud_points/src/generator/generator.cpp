#include "generator.h"

struct cloud_generator_data* new_cloud_generator (struct user_data *config)
{
    struct cloud_generator_data *generator = (struct cloud_generator_data*)malloc(sizeof(struct cloud_generator_data));

    generator->area = config->area;
    generator->num_points = config->num_points;

    generator->points = new std::vector<Point>();

    set_generator_function(generator,config);

    return generator;
}

void free_cloud_generator (struct cloud_generator_data *generator)
{
    delete generator->points;

    if (generator->handle)
        dlclose(generator->handle);

    free(generator);
}

void set_generator_function (struct cloud_generator_data *generator, struct user_data *config)
{
    std::string library_path(config->library_name->c_str());

    generator->handle = dlopen(library_path.c_str(),RTLD_LAZY);
    if (!generator->handle) 
    {
        fprintf(stderr,"%s\n",dlerror());
        exit(EXIT_FAILURE);
    }
    else
    {
        fprintf(stdout,"\n[generator] Cloud generator function library \"%s\" opened with sucess\n",library_path.c_str());
    }
    
    std::string function_name(config->function_name->c_str());

    generator->function = (set_cloud_generator_fn*)dlsym(generator->handle,function_name.c_str());
    if (dlerror() != NULL)  
    {
        fprintf(stderr, "[generator] \"%s\" function not found in the provided cloud generator function library\n",function_name.c_str());
        exit(EXIT_FAILURE);
    }
    else
    {
        fprintf(stdout, "[generator] Using the function \"%s\" as cloud generator function\n",function_name.c_str());
    }
    
}

void write_to_vtp (struct cloud_generator_data *generator)
{

    vtkSmartPointer<vtkAppendPolyData> appendFilter = vtkSmartPointer<vtkAppendPolyData>::New();
    
    // For each point create a sphere
    for (uint32_t i = 0; i < generator->points->size(); i++)
    {
        Point p = generator->points->at(i);

        vtkSmartPointer<vtkSphereSource> sphereSource = vtkSmartPointer<vtkSphereSource>::New();
        sphereSource->SetCenter(p.x,p.y,p.z);
        sphereSource->SetRadius(0.001*generator->area);
        sphereSource->Update();

        // Append the sphere to the filter
        appendFilter->AddInputConnection(sphereSource->GetOutputPort());
        appendFilter->Update();
    }

    // Write the file
    vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
    writer->SetFileName("output/cloud_points.vtp");
    writer->SetInputConnection(appendFilter->GetOutputPort());
    writer->Write();
}

void write_to_txt (struct cloud_generator_data *generator)
{
    std::vector<Point> *points = generator->points;

    FILE *file = fopen("output/cloud_points.txt","w+");

    fprintf(file,"%lu\n",points->size());
    for (uint32_t i = 0; i < points->size(); i++)    
    {
        Point p = generator->points->at(i);

        fprintf(file,"%g %g %g\n",p.x,p.y,p.z);
    }
        

    fclose(file);
}