#ifndef READER_H
#define READER_H

#include <cstdio>
#include <cstdlib>
#include <cstdbool>
#include <cstdint>
#include <cstring>

#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkUnstructuredGrid.h>
#include <vtkCell.h>
#include <vtkLine.h>
#include <vtkVertex.h>
#include <vtkHexahedron.h>
#include <vtkFloatArray.h>
#include <vtkCellData.h>
#include <vtkPointData.h>
#include <vtkGenericDataObjectReader.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkPolyDataWriter.h>

#include "../utils/utils.h"

class VTK_Reader
{
public:
    std::vector<Point> the_points;
public:
    VTK_Reader () { }
    VTK_Reader (std::string filename);
    uint32_t search_position (const double pos[]);
    vtkSmartPointer<vtkPolyData> convert_to_polydata ();
    bool check_duplicates ();
    void concatenate (VTK_Reader *input);
    void print ();
    void write (std::string filename);
};


void read_points (std::string filename, std::vector<Point> &the_points);

#endif