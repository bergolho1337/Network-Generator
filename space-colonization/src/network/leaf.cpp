#include "leaf.h"

Leaf::Leaf (const double x, const double y, const double z)
{
    this->is_reached = false;
    this->pos[0] = x;
    this->pos[1] = y;
    this->pos[2] = z;
}

void Leaf::print ()
{
    printf("(%g,%g,%g)\n",this->pos[0],this->pos[1],this->pos[2]);
}