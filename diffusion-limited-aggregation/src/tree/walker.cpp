#include "walker.h"

Walker::Walker ()
{
    sort_random_points(this->pos);
    this->stuck = false;
    this->radius = RADIUS;
}

Walker::Walker (const double x, const double y, const double z)
{
    this->pos[0] = x;
    this->pos[1] = y;
    this->pos[2] = z;

    this->stuck = false;
    this->radius = RADIUS;
}

void Walker::print ()
{
    printf("(%.5lf,%.5lf,%.5lf)\n",this->pos[0],this->pos[1],this->pos[2]);
}

// Generate a random point within the boundary of the domain
void sort_random_points (double pos[])
{
    uint8_t type = rand() % 4;
    double number = (double) rand() / (double) RAND_MAX;
    switch (type)
    {
        case 0: {
                    pos[0] = number * WIDTH;
                    pos[1] = 0.0;
                    pos[2] = 0.0;
                    break;
                }
        case 1: {
                    pos[0] = number * WIDTH;
                    pos[1] = HEIGHT;
                    pos[2] = 0.0;
                    break;
                }
        case 2: {
                    pos[0] = 0.0;
                    pos[1] = number * HEIGHT;
                    pos[2] = 0.0;
                    break;
                }
        case 3: {
                    pos[0] = WIDTH;
                    pos[1] = number * HEIGHT;
                    pos[2] = 0.0;
                    break;
                }
    }

}

void Walker::walk ()
{
    
    double dx = generate_random_number();
    double dy = generate_random_number();

    double new_x = this->pos[0] + dx;
    double new_y = this->pos[1] + dy;

    if (is_inside(new_x,new_y))
    {
        this->pos[0] += dx;
        this->pos[1] += dy;
    }
}

uint32_t Walker::is_stuck (std::vector<Walker*> the_tree)
{
    for (uint32_t i = 0; i < the_tree.size(); i++)
    {
        double d = calculate_distance(this->pos,the_tree[i]->pos);
        if (d < RADIUS * RADIUS * 4.0)
        {
            this->stuck = true;
            return i;
        }
    }
    return the_tree.size();
}

bool is_inside (const double x, const double y)
{
    if ( (x >= 0.0 && x <= WIDTH) && (y >= 0.0 && y <= HEIGHT) )
        return true;
    else
        return false;
}

double calculate_distance (const double p1[], const double p2[])
{
    return (p2[0] - p1[0])*(p2[0] - p1[0]) + (p2[1] - p1[1])*(p2[1] - p1[1]);
}

double generate_random_number ()
{
    uint8_t sign = rand() % 2;
    double value = (double) rand() / (double) RAND_MAX;
    if (sign)
        return value;
    else
        return -value;
}

void print_progress_bar (const uint32_t cur_iter, const uint32_t max_iter)
{
  double percentage = (cur_iter / (double) max_iter) * 100.0;
  uint32_t filled_length = nearbyint(100 * cur_iter / max_iter);
  
  std::string the_bar;
  for (uint32_t i = 0; i < filled_length; i++)
    the_bar += "\u2588";    // Using Unicode here ... (filled box)
  for (uint32_t i = 0; i < 100-filled_length; i++)
    the_bar += "-";
  
  printf("%s |%s| %.1lf%% %s\r","Progress",the_bar.c_str(),percentage,"Complete"); // Carriage return
  fflush(stdout); // Flush the standard output
}