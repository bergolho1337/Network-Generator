// Author: Lucas Berg

#include <iostream>
#include <vector>
#include <random>
#include <algorithm>

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <dirent.h>

using namespace std;

// ====================================================================================================================
// CONSTANTS AND MACROS
#define PRINT_LINE "--------------------------------------------------------------------------------------------"

void list_all_files_inside_folder (DIR *dir, const string folder_path)
{
    struct dirent *ent;

    printf("%s\n",PRINT_LINE);
    printf("Printing all the files inside the folder '%s'\n\n",folder_path.c_str());
    while ( (ent = readdir(dir)) )
    {
        printf("%s\n",ent->d_name);
    }
    printf("%s\n",PRINT_LINE);
}

void filter_files_inside_folder_by_extension (DIR *dir, const string folder_path, const string extension_name)
{
    struct dirent *ent;

    printf("%s\n",PRINT_LINE);
    printf("Printing files inside the folder '%s' with the extension '%s'\n\n",folder_path.c_str(),extension_name.c_str());
    while ( (ent = readdir(dir)) )
    {
        string filename = ent->d_name;
        if (filename.size() > 3)    // Avoid the '.' and '..'
        {
            string extension = filename.substr(filename.size()-3,filename.size()-1);

            if (extension_name == extension)
                printf("Filename = %s // Extension = %s\n",filename.c_str(),extension.c_str());

        }
    }
    printf("%s\n",PRINT_LINE);
}

void filter_and_store_files_inside_folder_by_extension (DIR *dir, const string folder_path, const string extension_name, vector<string> &files)
{
    struct dirent *ent;

    while ( (ent = readdir(dir)) )
    {
        string filename = ent->d_name;
        if (filename.size() > 3)    // Avoid the '.' and '..'
        {
            string extension = filename.substr(filename.size()-3,filename.size()-1);

            // Insert this file into the vector
            if (extension_name == extension)
                files.push_back(filename);

        }
    }
}

void read_files_numbers (vector<string> files, vector<uint32_t> &files_number)
{
    uint32_t start, end, offset;

    for (uint32_t i = 0; i < files.size(); i++)
    {
        start = files[i].find_last_of("_");
        end = files[i].find_first_of(".");
        offset = end - start - 1;

        string number = files[i].substr(start+1,offset);
        //printf("%s -- Start = %u // End = %u -- [Number = %s]\n",files[i].c_str(),start,end,number.c_str());
    
        size_t aux;
        files_number.push_back( stoi(number,&aux) );
    }
} 

void sort_files_by_number (vector<string> &files, vector<uint32_t> &files_number)
{
    for (uint32_t i = 0; i < files.size(); i++)
    {
        for (uint32_t j = 0; j < files.size(); j++)
        {
            if (files_number[j] > files_number[i])
            {
                uint32_t aux = files_number[i];
                files_number[i] = files_number[j];
                files_number[j] = aux;

                string aux2 = files[i];
                files[i] = files[j];
                files[j] = aux2;
            }
        }
    }
}

int main (int argc, char *argv[])
{
    if (argc-1 != 2)
    {
        printf("%s\n",PRINT_LINE);
        printf("Usage:> %s <input_folder> <extension_filter>\n",argv[0]);
        printf("%s\n",PRINT_LINE);
        printf("<input_folder> = Input folder with the files to be processed\n");
        printf("<extension_filter> = Name of the extension to be filter by the program\n");
        exit(EXIT_FAILURE);   
    }

    string folder_path = argv[1];
    string extension_name = argv[2];

    DIR *dir;

    dir = opendir(folder_path.c_str());
    if ( dir )
    {
        printf("[+] Folder '%s' has been open !\n",folder_path.c_str());

        // List all the files in folder
        //list_all_files_inside_folder(dir,folder_path);

        // Filter files with a specific extension
        //filter_files_inside_folder_by_extension(dir,folder_path,extension_name);

        // Filter files with a specific extension and sort them by their number
        vector<string> files;
        filter_and_store_files_inside_folder_by_extension(dir,folder_path,extension_name,files);

        vector<uint32_t> files_number;
        read_files_numbers(files,files_number);

        sort_files_by_number(files,files_number);

        for (uint32_t i = 0; i < files.size(); i++)
            printf("%s -- %u\n",files[i].c_str(),files_number[i]);

    }
    else
    {
        printf("[-] ERROR! Cannot open '%s' folder !\n",folder_path.c_str());
    }
    closedir(dir);

    return 0;
}
