#include "pmj.h"

PMJ::PMJ (const uint32_t index, const double x, const double y, const double z)
{
    this->index = index;
    this->pos[0] = x;
    this->pos[1] = y;
    this->pos[2] = z;
}

void PMJ::print ()
{
    printf("PMJ %u --> (%lf,%lf,%lf)\n",this->index,this->pos[0],this->pos[1],this->pos[2]);
}

void read_pmjs_from_file (vector<PMJ> &pmjs)
{
    uint32_t num_pmjs;
	double pos[3];

    char *pmj_filename = strdup("pmjs/private/elizabeth_LV_pmjs.pts");
	FILE *pmj_file = fopen(pmj_filename,"r");

    fscanf(pmj_file,"%u",&num_pmjs);
	for (uint32_t i = 0; i < num_pmjs; i++)
	{
		fscanf(pmj_file,"%lf %lf %lf",&pos[0],&pos[1],&pos[2]);

		PMJ pmj(i,pos[0],pos[1],pos[2]);
		pmjs.push_back(pmj);	
	}

	fclose(pmj_file);
}