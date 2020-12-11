// Author: Lucas Berg

#include <iostream>
#include <cstdio>
#include <cstdlib>

#include "cco/cco.h"

using namespace std;

int main (int argc, char *argv[])
{
    if (!(argc-1 == 1 || argc-1 == 4))
    {
        printf("-----------------------------------------------------------------------------------------------------------------------------------------------------------------------\n");
        printf("Usage:> %s <input_config> [additional_config]\n",argv[0]);
        printf("-----------------------------------------------------------------------------------------------------------------------------------------------------------------------\n");
        printf("Example:> %s input/test_LV.ini\n",argv[0]);
        printf("          %s inputs/test_RV_back_top.ini inputs/test_RV_front_top.ini inputs/test_RV_front_bottom.ini inputs/test_RV_back_bottom.ini\n",argv[0]);
        printf("-----------------------------------------------------------------------------------------------------------------------------------------------------------------------\n");
        exit(EXIT_FAILURE);
    }
    else
    {
        // Default case - Elizabeth LV
        if (argc-1 == 1)
        {
            User_Options *options = new User_Options(argv[1]);
            CCO_Network *the_network = new CCO_Network(options);

            the_network->grow_tree(options);

            delete the_network;
            delete options;
        }
        // Multiple networks case - Elizabeth RV
        else
        {
            std::vector<User_Options*> options_array;
            std::vector<CCO_Network*> network_array;

            // Define the terminals for the network linking
            double term_pos_1[3] = {54943.5,49193.5,12943.5};
            double term_pos_2[3] = {55943.5,48955.2,12693.5};
            double term_pos_3[3] = {55193.5,50193.5,17193.4};
            double term_pos_4[3] = {57757.5,48443.5,17193.4};
            double term_pos_5[3] = {57694,53943,27193};
            double term_pos_6[3] = {57943.5,50693.5,25883.5};

            //double link_1_offset = 0.556824;
            //double link_2_offset = 1.88936;
            //double link_3_offset = 1.55924;
            double link_1_offset = 0;
            double link_2_offset = 0;
            double link_3_offset = 0;

            for (uint32_t i = 0; i < 4; i++)
            {
                User_Options *options = new User_Options(argv[i+1]);
                options_array.push_back(options);

                CCO_Network *network = new CCO_Network(options);
                network_array.push_back(network);
            }

            // [LINKING]
            double his_lat_offset = network_array[0]->lat_offset;

            // Grow 'back_top' network and update the 'lat_offset' from the 'front_top' and 'front_bottom' networks
            network_array[0]->grow_tree(options_array[0]);
            Segment *term_1 = network_array[0]->get_terminal(term_pos_1);
            Segment *term_2 = network_array[0]->get_terminal(term_pos_3);
            network_array[1]->lat_offset = his_lat_offset + network_array[0]->calc_terminal_local_activation_time(term_1) + link_1_offset;
            network_array[2]->lat_offset = his_lat_offset + network_array[0]->calc_terminal_local_activation_time(term_2) + link_2_offset;

            // Grow 'front_top' and 'front_bottom' networks
            network_array[1]->grow_tree(options_array[1]);
            network_array[2]->grow_tree(options_array[2]);

            // Update the 'lat_offset' from the 'back_bottom'
            Segment *term_3 = network_array[2]->get_terminal(term_pos_5);
            network_array[3]->lat_offset = his_lat_offset + network_array[0]->calc_terminal_local_activation_time(term_2) + network_array[2]->calc_terminal_local_activation_time(term_3) + link_3_offset;

            // Grow 'back_bottom'
            network_array[3]->grow_tree(options_array[3]);

            // Get the reference to the other linking points
            Segment *term_4 = network_array[1]->get_terminal(term_pos_2);
            Segment *term_5 = network_array[2]->get_terminal(term_pos_4);
            Segment *term_6 = network_array[3]->get_terminal(term_pos_6);

            // Concatenate the networks
            CCO_Network *linked_network = network_array[0]->concatenate(network_array[1]);
            linked_network = linked_network->concatenate(network_array[2]);
            linked_network = linked_network->concatenate(network_array[3]);

            // Get the segments which have the linking points from the generated network
            term_1 = linked_network->get_terminal(term_pos_1);
            term_2 = linked_network->get_terminal(term_pos_2);
            term_3 = linked_network->get_terminal(term_pos_3);
            term_4 = linked_network->get_terminal(term_pos_4);
            term_5 = linked_network->get_terminal(term_pos_5);
            term_6 = linked_network->get_terminal(term_pos_6);

            // Link the segments
            //double dist_1 = euclidean_norm(term_1->dest->x,term_1->dest->y,term_1->dest->z,\
                                           term_2->src->x,term_2->src->y,term_2->src->z) * M_TO_UM;
            //double dist_2 = euclidean_norm(term_3->dest->x,term_3->dest->y,term_3->dest->z,\
                                           term_4->src->x,term_4->src->y,term_4->src->z) * M_TO_UM;
            //double dist_3 = euclidean_norm(term_5->dest->x,term_5->dest->y,term_5->dest->z,\
                                           term_6->src->x,term_6->src->y,term_6->src->z) * M_TO_UM;
            //printf("1--2 --> Dist=%g, LAT=%g\n",dist_1,dist_1/1900.0);
            //printf("3--4 --> Dist=%g, LAT=%g\n",dist_2,dist_2/1900.0);
            //printf("5--6 --> Dist=%g, LAT=%g\n",dist_3,dist_3/1900.0);
            printf("back_top --> LAT offset = %g\n",network_array[0]->lat_offset);
            printf("front_top --> LAT offset = %g\n",network_array[1]->lat_offset);
            printf("front_bottom --> LAT offset = %g\n",network_array[2]->lat_offset);
            printf("back_bottom --> LAT offset = %g\n",network_array[3]->lat_offset);
            linked_network->link_segments(term_1,term_2);
            linked_network->link_segments(term_3,term_4);
            linked_network->link_segments(term_5,term_6);

            // Set the 'output_dir' and adjust the radius
            linked_network->output_dir = "outputs/test_RV";
            linked_network->adjust_radius();
            
            printf("%s\n",PRINT_DOTS);
            printf("After linking ...\n");
            linked_network->print_network_info();
            linked_network->write_to_vtk_iteration();
            printf("%s\n",PRINT_DOTS);

            delete linked_network;

            for (uint32_t i = 0; i < 4; i++)
            {
                delete options_array[i];
                delete network_array[i];
            }
        }
    }
    

    return 0;
}
