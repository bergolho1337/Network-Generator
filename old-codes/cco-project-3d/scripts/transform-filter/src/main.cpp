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
        fprintf(stderr,"\t%s ../../saved-networks/elizabeth_minimize_volume_original_nterm600.vtk outputs/elizabeth_minimize_volume_original_nterm600_shifted.vtk\n",argv[0]);
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

        // ELIZABETH: Transform filter from CCO domain to MonoAlg3D domain !!!
    /*
        double t1[3] = {140,160,160};
        double t2[3] = {35,-40,40};
        double t3[3] = {20000,27500,27500};
        double t4[3] = {-1750,-3000,-3500};
        double t5[3] = {-1500,-2500,-3500};
        double s1[3] = {3750,3750,3750};
        double s2[3] = {0.25,0.25,0.25};
        double s3[3] = {1000,1000,1000};
        double s4[3] = {1.1,1.1,1.1};
        double s5[3] = {1.1,1.1,1.1};

        vtkSmartPointer<vtkPolyData> polydata_1 = transform_polydata(input_polydata,t1,180,s1);
        vtkSmartPointer<vtkPolyData> polydata_2 = transform_polydata(polydata_1,t2,180,s2);
        vtkSmartPointer<vtkPolyData> polydata_3 = transform_polydata(polydata_2,t3,0,s3);
        vtkSmartPointer<vtkPolyData> polydata_4 = transform_polydata(polydata_3,t4,0,s4);
        //vtkSmartPointer<vtkPolyData> polydata_5 = transform_polydata(polydata_4,t5,0,s5);
        
        vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New();
        writer->SetFileName(output_filename.c_str());
        writer->SetInputData(polydata_4);
        writer->Write();
    */

        // OXFORD
    /*
        double t1[3] = {0,0,0};
        double s1[3] = {1.0e+05,1.0e+05,1.0e+05};

        vtkSmartPointer<vtkPolyData> polydata_1 = transform_polydata(input_polydata,t1,0,s1);

        vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New();
        writer->SetFileName(output_filename.c_str());
        writer->SetInputData(polydata_1);
        writer->Write();    
    */
        // ELIZABETH BIVENTRICULAR
        double t1[3] = {-18226.2,-18226.2,-28226.2};
        double s1[3] = {4.0e+06,4.0e+06,4.0e+06};

        vtkSmartPointer<vtkPolyData> polydata_1 = transform_polydata(input_polydata,t1,0,s1);

        vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New();
        writer->SetFileName(output_filename.c_str());
        writer->SetInputData(polydata_1);
        writer->Write();    
    }
    else
    {
        fprintf(stderr,"[-] ERROR! The data should be Polydata!\n");
        exit(EXIT_FAILURE);
    }

  return EXIT_SUCCESS;
}
