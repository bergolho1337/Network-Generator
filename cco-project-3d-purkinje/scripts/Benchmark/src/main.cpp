/*
    Este programa recebe duas redes de Purkinje e calcula o erro nos tempos de ativacao nos PMJ's passados como entrada.
    Para isso, considera-se que a velocidade de propagacao na Purkinje eh constante.
*/

#include "reader/reader.h"
#include "graph/graph.h"

int main (int argc, char *argv[])
{
    if (argc-1 != 4)
    {
        printf("=============================================================================================================================\n");
        printf("Usage:> %s <input_ref_purkinje_network> <input_aprox_purkinje_network> <input_active_pmj> <output_filename>\n",argv[0]);
        printf("=============================================================================================================================\n");
        exit(EXIT_FAILURE);
    }

    std::string ref_filename = argv[1];
    std::string aprox_filename = argv[2];
    std::string pmj_filename = argv[3];
    std::string output_filename = argv[4];

    // [PURKINJE]
    VTK_Reader *ref_reader = new VTK_Reader(ref_filename);
    VTK_Reader *aprox_reader = new VTK_Reader(aprox_filename);

    Graph *ref_network = new Graph(ref_reader->the_points,ref_reader->the_lines);
    Graph *aprox_network = new Graph(aprox_reader->the_points,aprox_reader->the_lines);    

    ref_network->compute_errors(aprox_network,pmj_filename);
    //ref_network->write_terminals("outputs/reference_terminals.vtk");
    //ref_network->write_LAT("outputs/reference_LAT.vtk");
    
    aprox_network->write_LAT(output_filename.c_str());

    return 0;
}

/*
// Extra points for the RV cloud
    std::vector<uint32_t> ids;
    ids.push_back(791);
    ids.push_back(863);
    ids.push_back(862);
    ids.push_back(838);
    ids.push_back(887);
    ids.push_back(1133);
    ids.push_back(832);
    ids.push_back(782);
    ids.push_back(743);
    ids.push_back(1203);
    ids.push_back(1338);
    ids.push_back(1127);
    ids.push_back(1094);
    ids.push_back(1087);
    ids.push_back(1080);
    ids.push_back(1044);
    ids.push_back(1043);
    ids.push_back(467);
    ids.push_back(505);
    ids.push_back(616);
    ids.push_back(949);
    ids.push_back(976);
    ids.push_back(975);
    ids.push_back(948);
    ids.push_back(947);
    ids.push_back(1037);
    ids.push_back(2240);
    ids.push_back(2274);
    ids.push_back(2275);
    ids.push_back(1488);
    ids.push_back(1700);
    ids.push_back(1641);
    ids.push_back(1474);
    ids.push_back(1373);
    ids.push_back(323);
    ids.push_back(305);
    ids.push_back(279);
    ids.push_back(249);
    ids.push_back(333);
    ids.push_back(302);
    ids.push_back(205);
    ids.push_back(1844);
    ids.push_back(1985);
    ids.push_back(2038);
    ids.push_back(2159);
    ids.push_back(2259);
    ids.push_back(2208);
    ids.push_back(2203);
    ids.push_back(2256);
    ids.push_back(2197);
    ids.push_back(2195);
    ids.push_back(2033);
    ids.push_back(1796);
    ids.push_back(1512);
    ids.push_back(1665);
    ids.push_back(1873);
    ids.push_back(2126);
    ids.push_back(2253);
    ids.push_back(1273);
    ids.push_back(1346);
    ids.push_back(1384);
    ids.push_back(1383);
    ids.push_back(1403);
    ids.push_back(1663);
    ids.push_back(1662);
    ids.push_back(1450);
    ids.push_back(1512);
    ids.push_back(1665);
    ids.push_back(1873);
    ids.push_back(1659);

    printf("x y z\n");
    for (uint32_t i = 0; i < ids.size(); i++)
        printf("%g %g %g\n",ref_network->list_nodes[ids[i]].x,ref_network->list_nodes[ids[i]].y,ref_network->list_nodes[ids[i]].z);
    printf("%u\n",ids.size());
*/