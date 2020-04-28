#include <vtkSmartPointer.h>
#include <vtkArrowSource.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkTransform.h>
#include <vtkTransformFilter.h>
#include <vtkGenericDataObjectReader.h>
#include <vtkPolyDataWriter.h>
 
vtkSmartPointer<vtkPolyData> transform_polydata (vtkSmartPointer<vtkPolyData> polydata, const double t[], const double angle, const double s[])
{
    vtkSmartPointer<vtkTransform> transform = vtkSmartPointer<vtkTransform>::New();
    transform->Translate(t);
    transform->RotateY(angle);
    transform->Scale(s);

    vtkSmartPointer<vtkTransformFilter> transformFilter = vtkSmartPointer<vtkTransformFilter>::New();
    transformFilter->SetInputData(polydata);
    transformFilter->SetTransform(transform);

    vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->SetInputConnection(transformFilter->GetOutputPort());
    mapper->Update();

    vtkSmartPointer<vtkPolyData> output_polydata = vtkSmartPointer<vtkPolyData>::New();   
    output_polydata = mapper->GetInput();
    return output_polydata;
} 

int main (int argc, char *argv[])
{
    if (argc-1 != 2)
    {
        fprintf(stderr,"Usage:> %s <input_file> <output_file>\n",argv[0]);
        fprintf(stderr,"Example:\n");
        fprintf(stderr,"\t%s ../../saved-networks/elizabeth-purkinje-dla.vtk outputs/elizabeth-purkinje-dla_shifted.vtk\n",argv[0]);
        exit(EXIT_FAILURE);
    }

    std::string input_filename = argv[1];
    std::string output_filename = argv[2];

    // Get all data from the file
    vtkSmartPointer<vtkGenericDataObjectReader> reader = vtkSmartPointer<vtkGenericDataObjectReader>::New();
    reader->SetFileName(input_filename.c_str());
    reader->Update();

    // All of the standard data types can be checked and obtained like this:
    if(reader->IsFilePolyData())
    {
        vtkSmartPointer<vtkPolyData> input_polydata = reader->GetPolyDataOutput();

        // Transform filter from CCO domain to MonoAlg3D domain !!!
        double t2[3] = {35,-40,40};
        double t3[3] = {20000,27500,27500};
        double s2[3] = {0.25,0.25,0.25};
        double s3[3] = {1000,1000,1000};

        vtkSmartPointer<vtkPolyData> polydata_2 = transform_polydata(input_polydata,t2,180,s2);
        vtkSmartPointer<vtkPolyData> polydata_3 = transform_polydata(polydata_2,t3,0,s3);
        
        vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New();
        writer->SetFileName(output_filename.c_str());
        writer->SetInputData(polydata_3);
        writer->Write();
    }
    else
    {
        fprintf(stderr,"[-] ERROR! The data should be Polydata!\n");
        exit(EXIT_FAILURE);
    }

  return EXIT_SUCCESS;
}