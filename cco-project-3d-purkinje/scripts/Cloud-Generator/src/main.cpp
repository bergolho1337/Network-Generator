#include "reader/reader.h"

VTK_Reader* build_new_inactives_cloud (VTK_Reader *endocardium, VTK_Reader *actives, const double radius)
{
    VTK_Reader *result = new VTK_Reader();
    std::vector<bool> bitmask;
    bitmask.assign(endocardium->the_points.size(),false);
    for (uint32_t i = 0; i < actives->the_points.size(); i++)
    {
        //if (i == 0 || i == 4) continue;

        uint32_t counter = 0;
        double p1[3];
        p1[0] = actives->the_points[i].x;
        p1[1] = actives->the_points[i].y;
        p1[2] = actives->the_points[i].z;

        for (uint32_t j = 0; j < endocardium->the_points.size(); j++)
        {
            double p2[3];
            p2[0] = endocardium->the_points[j].x;
            p2[1] = endocardium->the_points[j].y;
            p2[2] = endocardium->the_points[j].z;
            
            double dist = sqrt(powf(p1[0]-p2[0],2)+powf(p1[1]-p2[1],2)+powf(p1[2]-p2[2],2));
            if (dist < radius && !bitmask[j])
            {
                if (counter % 5 == 0) result->the_points.push_back(endocardium->the_points[j]);
                bitmask[j] = true;
                counter++;
            }
        }
    }

    // Check for duplicates
    if (result->check_duplicates()) exit(EXIT_FAILURE);

    return result;
}

int main (int argc, char *argv[])
{
    if (argc-1 != 4)
    {
        printf("=============================================================================================================================\n");
        printf("Usage:> %s <input_inactives_cloud> <input_actives_cloud> <input_endocardium_cloud> <output_inactives_cloud>\n",argv[0]);
        printf("=============================================================================================================================\n");
        exit(EXIT_FAILURE);
    }

    std::string inactive_filename = argv[1];
    std::string active_filename = argv[2];
    std::string endocardium_filename = argv[3];
    std::string output_filename = argv[4];

    VTK_Reader *inactive_cloud = new VTK_Reader(inactive_filename);
    VTK_Reader *active_cloud = new VTK_Reader(active_filename);
    VTK_Reader *endocardium_cloud = new VTK_Reader(endocardium_filename);

    VTK_Reader *output_cloud = build_new_inactives_cloud(endocardium_cloud,active_cloud,0.001);
    inactive_cloud->concatenate(output_cloud);
    inactive_cloud->write(output_filename);
    
    delete inactive_cloud;
    delete active_cloud;
    delete endocardium_cloud;

    return 0;
}
