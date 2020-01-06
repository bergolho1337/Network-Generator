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

                struct point p;
                p.x = pos[0];
                p.y = pos[1];
                p.z = pos[2];
                //printf("Bifurcation position test -> (%g,%g,%g)\n",pos[0],pos[1],pos[2]);
                
                test_positions.push_back(p);

            }
        }
    }
}

SET_LOCAL_OPTIMIZATION_FUNCTION (rafael_local_optimization)
{

    const double d_eta = 1.0 / NE;
    const double d_neta = 1.0 / NE;
    const uint32_t npts = NE + 1;

    // Get a reference to the 4 vertices from the bifurcation
    struct point *A = ibiff->value->src->value;
    struct point *B = iconn->value->dest->value;
    struct point *Bif = ibiff->value->dest->value;
    struct point *T = inew->value->dest->value;

    // Calculate the parametrization
    double t = 0.3; // ???
    double E[3], F[3], G[3];
    
    E[0] = B->x + (A->x - B->x) * (1-t);
    E[1] = B->y + (A->y - B->y) * (1-t);
    E[2] = B->z + (A->z - B->z) * (1-t);

    F[0] = T->x + (Bif->x - T->x) * (t);
    F[1] = T->y + (Bif->y - T->y) * (t);
    F[2] = T->z + (Bif->z - T->z) * (t);

    G[0] = B->x + (A->x - B->x) * (t);
    G[1] = B->y + (A->y - B->y) * (t);
    G[2] = B->z + (A->z - B->z) * (t);
    
    // DEBUG
    //printf("G1 = (%g,%g,%g)\n",G[0],G[1],G[2]);
    //printf("G2 = (%g,%g,%g)\n",F[0],F[1],F[2]);
    //printf("G3 = (%g,%g,%g)\n",E[0],E[1],E[2]);

    // Build the test points for the local optimization based on the triangle composed by
    // 'G1', 'G2' and 'G3'
    for (uint32_t i = npts; i >= 1; i--)
    {
        for (uint32_t j = 1; j <= i; j++)
        {
            double neta = (npts - i) * d_neta;
            double eta = (j - 1) * d_eta;
            double phi = (1 - eta - neta);

            double pos[3];
            pos[0] = (phi * G[0]) + (eta * F[0]) + (neta * E[0]);
            pos[1] = (phi * G[1]) + (eta * F[1]) + (neta * E[1]);
            pos[2] = (phi * G[2]) + (eta * F[2]) + (neta * E[2]);

            printf("Bifurcation position test -> (%g,%g,%g)\n",pos[0],pos[1],pos[2]);

            struct point p;
            p.x = pos[0];
            p.y = pos[1];
            p.z = pos[2];
            
            test_positions.push_back(p);
        }
    }
}
