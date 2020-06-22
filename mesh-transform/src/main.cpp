// ------------------------------------------------------------------
// Example that demonstrates the use of the 'getopt.h' library
// It is used to parse command line options
// ------------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>     // Library used to parse command line options

#include "vtk_utils/vtk_utils.h"

#define PRINT_LINE "======================================================================"

// List of command line options. If the letter option is follow by a colon them the option requires a an argument.
static const char *opt_string =   "i:o:r:t:s:gh";

static const struct option long_options[] = {
        { "input", required_argument, NULL, 'i'},
        { "output", required_argument, NULL, 'o'},
        { "rotate", required_argument, NULL, 'r'},
        { "translate", required_argument, NULL, 't'},
        { "scale", required_argument, NULL, 's'},
        { "origin", no_argument, NULL, 'g'},
        { "help", no_argument, NULL, 'h' },
        { NULL, no_argument, NULL, 0 }
};

void usage (const char pName[]) {

    printf("%s\n",PRINT_LINE);
    printf("Usage:> %s -option\n",pName);
    printf("option =\n");
    printf("\t-i <name_of_file> = input file\n");
    printf("\t-o <name_of_file> = output file\n");
    printf("\t-r \"x y z\" = rotate mesh {in degrees}\n");
    printf("\t-t \"dx dy dz\" = translate mesh\n");
    printf("\t-s \"x y z\" = scale mesh\n");
    printf("\t-g  = translate mesh to origin {0,0,0}\n");
    printf("\t-h = help\n");
    printf("%s\n",PRINT_LINE);
    exit(EXIT_FAILURE);
    
}

int main (int argc, char *argv[]) {

    int opt;
    int option_index = 0;

    // Get the command line options based on the configured options:
    opt = getopt_long(argc,argv,opt_string,long_options,&option_index);
    if (opt == -1)
        usage(argv[0]);

    std::string input_filename;
    bool have_input_file = false;
    std::string output_filename;
    bool have_output_file = false;
    std::string rotate_config;
    bool have_to_rotate = false;
    std::string translate_config;
    bool have_to_translate = false;
    bool have_to_translate_to_origin = false;
    std::string scale_config;
    bool have_to_scale = false;

    while (opt != -1) {

        switch (opt) {

            case 'i': input_filename = optarg;
                      have_input_file = true;
                      printf("[+] Input file '%s' ...\n",optarg);
                      break;
            case 'o': output_filename = optarg;
                      have_output_file = true;
                      printf("[+] Output file '%s' ...\n",optarg);
                      break;
            case 'r': rotate_config = optarg;
                      have_to_rotate = true;
                      printf("[+] Rotating mesh by (%s) degrees ...\n",optarg);
                      break;
            case 't': translate_config = optarg;
                      have_to_translate = true;
                      printf("[+] Translating mesh by (%s) ...\n",optarg);
                      break;
            case 's': scale_config = optarg;
                      have_to_scale = true;
                      printf("[+] Scaling mesh by (%s) ...\n",optarg);
                      break;
            case 'g': have_to_translate_to_origin = true;
                      printf("[+] Translating mesh to origin ...\n");
                      break;
            case 'h': usage(argv[0]);
                      break;
            default:  printf("[-] ERROR ! Invalid option\n");
                      break;

        }
        opt = getopt_long(argc,argv,opt_string,long_options,&option_index);
    }

    if (have_input_file) {

        std::string extension = get_extension(input_filename);

        // Supress warnings ...
        vtkObject::GlobalWarningDisplayOff();

        vtkSmartPointer<vtkPolyData> polydata_grid = vtkSmartPointer<vtkPolyData>::New();
        vtkSmartPointer<vtkUnstructuredGrid> unstructured_grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
        bool is_polydata = false;
        bool is_unstructured = false;

        if (extension == "vtu") 
        {
            vtkSmartPointer<vtkXMLUnstructuredGridReader> reader = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
            reader->SetFileName(input_filename.c_str());
            reader->Update();
            
            unstructured_grid = reader->GetOutput();
            is_unstructured = true;
            printf("Number of points = %u\n",unstructured_grid->GetNumberOfPoints());
        }
        else if (extension == "vtp") 
        {
            vtkSmartPointer<vtkXMLPolyDataReader> reader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
            reader->SetFileName(input_filename.c_str());
            reader->Update();
            
            polydata_grid = reader->GetOutput();
            is_polydata = true;
            printf("Number of points = %u\n",polydata_grid->GetNumberOfPoints());
        }
        else if (extension == "vtk") 
        {
            vtkSmartPointer<vtkGenericDataObjectReader> reader = vtkSmartPointer<vtkGenericDataObjectReader>::New();
            reader->SetFileName(input_filename.c_str());
            reader->Update();

            if(reader->IsFilePolyData())
            {
                polydata_grid = reader->GetPolyDataOutput();
                is_polydata = true;
                printf("Number of points = %u\n",polydata_grid->GetNumberOfPoints());

/*
                for (int i = 0; i < polydata_grid->GetNumberOfCells(); i++)
                {
                    vtkCell *cell = polydata_grid->GetCell(i);
                    printf("Cell %d -- Type = %d\n",i,cell->GetCellType());
                    
                    //vtkLine *line = dynamic_cast<vtkLine*>(cell);

                    //int src = line->GetPointIds()->GetId(0);
                    //int dest = line->GetPointIds()->GetId(1);

                    //cout << "\t " << src << " --> " << dest << endl;
                }
*/
                
            }
            else if (reader->IsFileUnstructuredGrid())
            {
                unstructured_grid = reader->GetUnstructuredGridOutput();
                is_unstructured = true;
                printf("Number of points = %u\n",unstructured_grid->GetNumberOfPoints());
            }

        }
        else 
        {
            fprintf(stderr,"[-] ERROR! Invalid input file extension '%s'! You must provide an vtk/vtu/vtp file format!\n",extension.c_str());
            exit(EXIT_FAILURE);
        }

        vtkSmartPointer<vtkPolyData> output_polydata = vtkSmartPointer<vtkPolyData>::New();
        output_polydata = polydata_grid;

        if (have_to_translate) {

            double disp[3];
            std::stringstream ss(translate_config);
            ss >> disp[0] >> disp[1] >> disp[2];

            if (is_polydata)
                output_polydata = translate_polydata(output_polydata,disp);
            else if (is_unstructured)
                printf("[-] Need to implement!\n");

        }

        if (have_to_rotate) {
            
            double rot[3];
            std::stringstream ss(rotate_config);
            ss >> rot[0] >> rot[1] >> rot[2];

            output_polydata = rotate_polydata(output_polydata,rot);
        }

        if (have_to_scale) {
            
            double scale[3];
            std::stringstream ss(scale_config);
            ss >> scale[0] >> scale[1] >> scale[2];

            output_polydata = scale_polydata(output_polydata,scale);
        }

        if (have_to_translate_to_origin) {

            double bounds[6];
            output_polydata->GetBounds(bounds);
            
            double min[3];
            min[0] = -bounds[0];
            min[1] = -bounds[2];
            min[2] = -bounds[4];
            printf("[+] Translating (%g,%g,%g)\n",min[0],min[1],min[2]);

            output_polydata = translate_polydata(output_polydata,min);
        }

        if (have_output_file) {

            if (is_polydata) {

                vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
                writer->SetFileName(output_filename.c_str());
                writer->SetInputData(output_polydata);
                writer->Write();
            }
            else if (is_unstructured) {

                printf("[-] Need to implement!\n");
            }
            
            
        }
        else {
            fprintf(stderr,"[-] ERROR! You must provide an output file!\n");
            exit(EXIT_FAILURE);
        }

        
    }
    else {
        fprintf(stderr,"[-] ERROR! You must provide an input file!\n");
        exit(EXIT_FAILURE);
    }

    

    return 0;
}