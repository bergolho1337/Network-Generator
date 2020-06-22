#ifndef VTK_UTILS_H
#define VTK_UTILS_H

#include <iostream>
#include <string>
#include <sstream>

#include <vtkGenericDataObjectReader.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkCylinderSource.h>
#include <vtkStructuredGrid.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkTransform.h>
#include <vtkTransformFilter.h>
#include <vtkPolyDataMapper.h>
#include <vtkPolyDataWriter.h>
#include <vtkAppendFilter.h>
#include <vtkDataSetMapper.h>

vtkSmartPointer<vtkPolyData> translate_polydata (vtkPolyData *polydata, const double disp[]);
vtkSmartPointer<vtkPolyData> rotate_polydata (vtkPolyData *polydata, const double rot[]);
vtkSmartPointer<vtkPolyData> scale_polydata (vtkPolyData *polydata, const double scale[]);

std::string get_extension (std::string filename);

#endif