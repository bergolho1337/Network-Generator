#include "vtk_utils.h"

std::string get_extension (std::string filename) {

    size_t len = filename.size();
    size_t pos = filename.find(".");

    std::string extension = filename.substr(pos+1,len);

    return extension;
}

vtkSmartPointer<vtkPolyData> translate_polydata (vtkPolyData *polydata, const double disp[]) {

    vtkSmartPointer<vtkTransform> transform = vtkSmartPointer<vtkTransform>::New();
    transform->Translate(disp);

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

vtkSmartPointer<vtkPolyData> rotate_polydata (vtkPolyData *polydata, const double rot[]) {

    vtkSmartPointer<vtkTransform> transform = vtkSmartPointer<vtkTransform>::New();
    transform->RotateX(rot[0]);
    transform->RotateY(rot[1]);
    transform->RotateZ(rot[2]);

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

vtkSmartPointer<vtkPolyData> scale_polydata (vtkPolyData *polydata, const double scale[]) {

    vtkSmartPointer<vtkTransform> transform = vtkSmartPointer<vtkTransform>::New();
    transform->Scale(scale);

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
