#include "reader.h"

uint32_t counter_points_table = 0;

void insert_point_into_table (map<Point_Table,uint32_t> &point_table, const double pos[])
{
    map<Point_Table,uint32_t>::iterator it;
    Point_Table point(pos[0],pos[1],pos[2]);

    it = point_table.find(point);
    if (it == point_table.end())
    {
        point_table.insert(pair<Point_Table,uint32_t>(point,counter_points_table));
        counter_points_table++;
    }
}

void read_faces (const char *filename, vector<Point> &points, vector<Face> &faces)
{
    vtkSmartPointer<vtkSTLReader> reader = vtkSmartPointer<vtkSTLReader>::New();
    reader->SetFileName(filename);
    reader->Update();

    // Visualize
    vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->SetInputConnection(reader->GetOutputPort());
    mapper->Update();

        // Get the reference to the Polydata
    vtkPolyData *polydata_grid = mapper->GetInput();

    // Read cells
    int num_cells = polydata_grid->GetNumberOfCells();
    int num_points = polydata_grid->GetNumberOfPoints();

    //cout << "Number of points = " << num_points << endl;
    //cout << "Number of cells = " << num_cells << endl;

    // Read the points
    for (int i = 0; i < num_points; i++)
    {
        double pos[3];

        polydata_grid->GetPoint(i,pos);

        Point point(i,pos[0],pos[1],pos[2]);

        points.push_back(point);
        //cout << "Point " << i << " = (" << pos[0] << "," << pos[1] << "," << pos[2] << ")" << endl;
    }

    // Read the cells
    for (int i = 0; i < num_cells; i++)
    {
        //cout << "Cell " << i << endl;

        vtkCell *cell = polydata_grid->GetCell(i);

        vtkTriangle *line = dynamic_cast<vtkTriangle*>(cell);

        uint32_t vertex_1 = line->GetPointIds()->GetId(0);
        uint32_t vertex_2 = line->GetPointIds()->GetId(1);
        uint32_t vertex_3 = line->GetPointIds()->GetId(2);

        Face face(vertex_1,vertex_2,vertex_3);
        
        faces.push_back(face);
        //cout << "\t " << vertex_1 << " " << vertex_2 << " " << vertex_3 << endl;
    }
}

void read_points_from_pts (const char *filename, vector<Point> &points)
{
    FILE *file = fopen(filename,"r");

    uint32_t num_points;
    fscanf(file,"%u",&num_points);

    for (uint32_t i = 0; i < num_points; i++)
    {
        double x, y, z;
        fscanf(file,"%lf %lf %lf",&x,&y,&z);

        Point p(i,x,y,z);
        points.push_back(p);
    }

    fclose(file);
}

void print_points (vector<Point> points)
{
    for (int i = 0; i < (int)points.size(); i++)
        points[i].print();
}

void print_faces (vector<Face> faces)
{
    for (int i = 0; i < (int)faces.size(); i++)
    {
        printf("[Face %d]\n",i);
        faces[i].print();
    }
}

void write_marked_faces_to_stl(vector<Point> points, vector<Face> faces)
{

    vtkSmartPointer<vtkPoints> the_points = vtkSmartPointer<vtkPoints>::New();
    for (uint32_t i = 0; i < points.size(); i++)
        the_points->InsertNextPoint(points[i].x,points[i].y,points[i].z);
    
    vtkSmartPointer<vtkCellArray> the_cell_array = vtkSmartPointer<vtkCellArray>::New();
    for (uint32_t i = 0; i < faces.size(); i++)
    {
        uint32_t vertex_index_1 = faces[i].vertex_index_1;
        uint32_t vertex_index_2 = faces[i].vertex_index_2;
        uint32_t vertex_index_3 = faces[i].vertex_index_3;

        if (faces[i].marked)
        {
            vtkSmartPointer<vtkTriangle> triangle = vtkSmartPointer<vtkTriangle>::New();

            triangle->GetPointIds()->SetId(0,vertex_index_1);
            triangle->GetPointIds()->SetId(1,vertex_index_2);
            triangle->GetPointIds()->SetId(2,vertex_index_3);

            the_cell_array->InsertNextCell(triangle);
        }
    }

    vtkSmartPointer<vtkPolyData> polydata_grid = vtkSmartPointer<vtkPolyData>::New();
    polydata_grid->SetPoints(the_points);
    polydata_grid->SetPolys(the_cell_array);
    
    vtkSmartPointer<vtkSTLWriter> stlWriter = vtkSmartPointer<vtkSTLWriter>::New();
    stlWriter->SetFileName("outputs/processed_mesh.stl");
    stlWriter->SetInputData(polydata_grid);
    stlWriter->Write();

/*
    FILE *file = fopen("outputs/processed_mesh.stl","w+");

    fprintf(file,"solid ascii\n");
    for (uint32_t i = 0; i < faces.size(); i++)
    {
        if (faces[i].marked)
        {
            write_face_to_stl(file,faces[i],points);
        }
    }
    fprintf(file,"endsolid\n");

    fclose(file);
*/
}

void write_face_to_stl(FILE *file, Face face, vector<Point> points)
{
    uint32_t vertex_index_1 = face.vertex_index_1;
    uint32_t vertex_index_2 = face.vertex_index_2;
    uint32_t vertex_index_3 = face.vertex_index_3;
    
    // TODO: Maybe include the original normal here ...
    fprintf(file," facet normal 0 0 1\n");
    fprintf(file,"  outer loop\n");
    fprintf(file,"   vertex %g %g %g\n",points[vertex_index_1].x,points[vertex_index_1].y,points[vertex_index_1].z);
    fprintf(file,"   vertex %g %g %g\n",points[vertex_index_2].x,points[vertex_index_2].y,points[vertex_index_2].z);
    fprintf(file,"   vertex %g %g %g\n",points[vertex_index_3].x,points[vertex_index_3].y,points[vertex_index_3].z);
    fprintf(file,"  endloop\n");
    fprintf(file," endfacet\n");

}