#include "local_optimization.h"

SET_LOCAL_OPTIMIZATION_FUNCTION (default_local_optimization)
{

    double delta_e = 1.0 / NE;

    // Get a reference to the 3 vertices surrounding the triangle area
    struct point *g1 = ibiff->value->src->value;
    struct point *g2 = iconn->value->dest->value;
    struct point *g3 = inew->value->dest->value;

    //printf("G1 = (%g,%g,%g)\n",g1->x,g1->y,g1->z);
    //printf("G2 = (%g,%g,%g)\n",g2->x,g2->y,g2->z);
    //printf("G3 = (%g,%g,%g)\n",g3->x,g3->y,g3->z);

    // Build the test points for the local optimization based on the triangle composde by
    // 'iconn', 'ibiff' and 'inew'
    for (uint32_t i = 0; i <= NE; i++)
    {
        for (uint32_t j = 0; j <= NE-i; j++)
        {
            double epsilon = i*delta_e;
            double eta = j*delta_e;

            if (!is_corner(i,j,NE))
            {
                // Build the phi array
                double phi[3];
                phi[0] = 1.0 - epsilon - eta;
                phi[1] = epsilon;
                phi[2] = eta;

                double pos[3];
                pos[0] = (phi[0]*g1->x) + (phi[1]*g2->x) + (phi[2]*g3->x);
                pos[1] = (phi[0]*g1->y) + (phi[1]*g2->y) + (phi[2]*g3->y);
                pos[2] = (phi[0]*g1->z) + (phi[1]*g2->z) + (phi[2]*g3->z);

                struct point *p = new_point(pos);
                test_positions.push_back(p);

            }
        }
    }
}
