#include "leaf.h"

Leaf::Leaf (const uint32_t width, const uint32_t height)
{
    this->is_reached = false;
    this->pos[0] = generate_random_number() * width;
    this->pos[1] = generate_random_number() * height - 100;
    this->pos[2] = 0.0;
}

void Leaf::print ()
{
    printf("(%g,%g,%g)\n",this->pos[0],this->pos[1],this->pos[2]);
}