#include "generator.h"

Cloud_Generator::Cloud_Generator (User_Data *config)
{
    this->num_points = config->num_points;

    this->points = new std::vector<Point>();

    set_generator_function(config);
}

Cloud_Generator::~Cloud_Generator ()
{
    delete this->points;

    if (this->handle)
        dlclose(this->handle);
}

void Cloud_Generator::set_generator_function (User_Data *config)
{
    std::string library_path(config->library_name->c_str());

    this->handle = dlopen(library_path.c_str(),RTLD_LAZY);
    if (!this->handle) 
    {
        fprintf(stderr,"%s\n",dlerror());
        exit(EXIT_FAILURE);
    }
    else
    {
        fprintf(stdout,"\n[generator] Cloud generator function library \"%s\" opened with sucess\n",library_path.c_str());
    }
    
    std::string function_name(config->function_name->c_str());

    this->function = (set_cloud_generator_fn*)dlsym(this->handle,function_name.c_str());
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

void Cloud_Generator::write_to_pts ()
{
    fprintf(stdout,"[generator] Writing cloud of points in PTS format ...\n");

    std::vector<Point> *points = this->points;

    FILE *file = fopen("output/cloud_points.pts","w+");

    fprintf(file,"%lu\n",points->size());
    for (uint32_t i = 0; i < points->size(); i++)    
    {
        Point p = this->points->at(i);

        fprintf(file,"%g %g %g\n",p.x,p.y,p.z);
    }
        
    fclose(file);
}

void Cloud_Generator::write_to_vtk ()
{
    fprintf(stdout,"[generator] Writing cloud of points in VTK format ...\n");

    std::vector<Point> *points = this->points;

    FILE *file = fopen("output/cloud_points.vtk","w+");

    fprintf(file,"# vtk DataFile Version 3.0\n");
    fprintf(file,"Cloud\n");
    fprintf(file,"ASCII\n");
    fprintf(file,"DATASET POLYDATA\n");
    fprintf(file,"POINTS %lu float\n",points->size());

    for (uint32_t i = 0; i < points->size(); i++)    
    {
        Point p = this->points->at(i);

        fprintf(file,"%g %g %g\n",p.x,p.y,p.z);
    }

    fprintf(file,"VERTICES %lu %lu\n",points->size(),points->size()*2);    
    for (uint32_t i = 0; i < points->size(); i++)
        fprintf(file,"1 %u\n",i);

    fclose(file);
}

void Cloud_Generator::write_to_vtp ()
{
    fprintf(stdout,"[generator] Writing cloud of points in VTP format ...\n");

    vtkSmartPointer<vtkAppendPolyData> appendFilter = vtkSmartPointer<vtkAppendPolyData>::New();
    
    // For each point create a sphere
    for (uint32_t i = 0; i < this->points->size(); i++)
    {
        Point p = this->points->at(i);

        vtkSmartPointer<vtkSphereSource> sphereSource = vtkSmartPointer<vtkSphereSource>::New();
        sphereSource->SetCenter(p.x,p.y,p.z);
        //sphereSource->SetRadius(0.001*generator->area);
	    sphereSource->SetRadius(1.0);        
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

