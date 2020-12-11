#include "local_optimization.h"

#include "../cco/cco.h"
#include "../options/local_optimization_config.h"

void LocalOptimization::move_bifurcation_location (Segment *iconn, Segment *ibiff, Segment *inew, Point *p)
{
    ibiff->dest->x = p->x;
    ibiff->dest->y = p->y;
    ibiff->dest->z = p->z;

    iconn->src->x = p->x;
    iconn->src->y = p->y;
    iconn->src->z = p->z;
    
    inew->src->x = p->x;
    inew->src->y = p->y;
    inew->src->z = p->z;
}

void LocalOptimization::move_bifurcation_location (Segment *iconn, Segment *ibiff, Segment *inew, const double pos[])
{
    ibiff->dest->x = pos[0];
    ibiff->dest->y = pos[1];
    ibiff->dest->z = pos[2];

    iconn->src->x = pos[0];
    iconn->src->y = pos[1];
    iconn->src->z = pos[2];
    
    inew->src->x = pos[0];
    inew->src->y = pos[1];
    inew->src->z = pos[2];
}

void LocalOptimization::free_test_positions (std::vector<Point*> &test_positions)
{
    for (uint32_t i = 0; i < test_positions.size(); i++)
    {
        delete test_positions[i];
        test_positions[i] = NULL;
    }
}

void DefaultOptimization::optimize (Segment *iconn, Segment *ibiff, Segment *inew,\
                                    std::vector<Point*> &test_positions)
{
    // implement me ...
}

void RafaelOptimization::optimize (Segment *iconn, Segment *ibiff, Segment *inew,\
                                    std::vector<Point*> &test_positions)
{
    const double d_eta = 1.0 / NE;
    const double d_neta = 1.0 / NE;
    const uint32_t npts = NE + 1;

    // Get a reference to the 4 vertices from the bifurcation
    Point *A = ibiff->src;
    Point *B = iconn->dest;
    Point *Bif = ibiff->dest;
    Point *T = inew->dest;

    // Calculate the parametrization
    double t = 0.3; // Why this value ???
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
    
    // Build the test points for the local optimization based on the triangle composed by
    // 'G1', 'G2' and 'G3'
    uint32_t counter = 0;
    for (uint32_t i = npts; i >= 1; i--)
    {
        for (uint32_t j = 1; j <= i; j++)
        {
            if (!is_vertex(i,j))
            {
                double neta = (npts - i) * d_neta;
                double eta = (j - 1) * d_eta;
                double phi = (1 - eta - neta);

                double pos[3];
                pos[0] = (phi * G[0]) + (eta * F[0]) + (neta * E[0]);
                pos[1] = (phi * G[1]) + (eta * F[1]) + (neta * E[1]);
                pos[2] = (phi * G[2]) + (eta * F[2]) + (neta * E[2]);

                Point *p = new Point(counter);
                p->setCoordinate(pos);
                p->setActive(false);
                p->setLAT(0.0);
                
                test_positions.push_back(p);
                counter++;
            }
        }
    }
}

bool RafaelOptimization::is_vertex (const uint32_t i, const uint32_t j)
{
    return ((i == (NE+1) && j == (NE+1)) || (i == (NE+1) && j == 1) || (i == 1 && j == 1)) ? true : false;
}