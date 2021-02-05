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
#include <vtkHexahedron.h>
#include <vtkVertex.h>
#include <vtkFloatArray.h>
#include <vtkCellData.h>
#include <vtkPointData.h>
#include <vtkFieldData.h>
#include <vtkGenericDataObjectReader.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkDijkstraGraphGeodesicPath.h>
#include <vtkSTLReader.h>

#include "../utils/utils.h"

class VTK_Reader
{
public:
    std::vector<Point> the_points;
public:
    VTK_Reader () { }
    VTK_Reader (std::string filename);
    void concatenate (VTK_Reader *input);
    void print ();
    void write (std::string filename);
};

void read_terminal_positions (std::string filename, std::vector<Point> &terminals);

#endif