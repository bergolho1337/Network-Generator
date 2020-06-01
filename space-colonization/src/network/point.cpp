#include "point.h"

Point::Point (const uint32_t index, const double x, const double y, const double z)
{
    this->index = index;
    this->pos[0] = x;
    this->pos[1] = y;
    this->pos[2] = z;
}

void Point::print ()
{
    printf("(%g,%g,%g)\n",this->pos[0],this->pos[1],this->pos[2]);
}