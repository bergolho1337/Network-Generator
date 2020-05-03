#include "branch.h"

Branch::Branch (const uint32_t id, const double x, const double y, const double z,\
            const double dx, const double dy, const double dz,\
            uint32_t parent)
{
    this->id = id;

    this->pos[0] = x;
    this->pos[1] = y;
    this->pos[2] = z;
    
    this->dir[0] = dx;
    this->dir[1] = dy;
    this->dir[2] = dz;
    
    this->original_dir[0] = dx;
    this->original_dir[1] = dy;
    this->original_dir[2] = dz;

    this->length = 5.0;
    this->counter = 0;
    this->parent = parent;
}

void Branch::get_next_branch_position (double new_pos[])
{
    for (uint32_t i = 0; i < 3; i++)
        new_pos[i] = this->pos[i] + this->dir[i]*this->length;
}

void Branch::reset ()
{
    this->counter = 0;
    for (uint32_t i = 0; i < 3; i++)
        this->dir[i] = this->original_dir[i];
}

void Branch::print ()
{
    printf("Position = (%g,%g,%g) -- Direction = (%g,%g,%g)\n",this->pos[0],this->pos[1],this->pos[2],this->dir[0],this->dir[1],this->dir[2]);
}