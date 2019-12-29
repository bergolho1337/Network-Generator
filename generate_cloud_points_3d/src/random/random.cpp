#include "random.h"

Random_Generator::Random_Generator (const uint32_t seed)
{
    this->counter = 0;
    this->array = NULL;
    this->seed = seed;
}

Random_Generator::~Random_Generator ()
{
    if (this->array)
        delete [] this->array;
}

void Random_Generator::generate ()
{
    this->array = new double[RANDOM_ARRAY_SIZE];

    dsfmt_t dsfmt;
    dsfmt_init_gen_rand(&dsfmt, this->seed);
    for (uint32_t i = 0; i < RANDOM_ARRAY_SIZE; ++i) 
    {
        this->array[i] = dsfmt_genrand_close_open(&dsfmt);
    }
}

void Random_Generator::print ()
{
    for (uint32_t i = 0; i < RANDOM_ARRAY_SIZE; ++i)
    {
        printf("%g\n",this->array[i]);
    }
}

double Random_Generator::get_value ()
{
    double value = this->array[this->counter];

    this->counter++;
    if (this->counter > RANDOM_ARRAY_SIZE)
        this->counter = 0;

    return value;
}