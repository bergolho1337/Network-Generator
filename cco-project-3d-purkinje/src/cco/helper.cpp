#include "helper.h"

bool check_bifurcation_rule (FILE *log_file, const double gamma, std::vector<Segment*> s_list)
{
    for (uint32_t i = 0; i < s_list.size(); i++)
    {
        Segment *cur_segment = s_list[i];
        Segment *left = cur_segment->left;
        Segment *right = cur_segment->right;

        if (left != NULL && right != NULL)
        {
            double r = pow(cur_segment->radius,gamma);
            double r_left = pow(left->radius,gamma);
            double r_right = pow(right->radius,gamma);

            double diff = r_left + r_right - r;
            fprintf(log_file,"%g = %g + %g --> %g = %g\n",r,r_left,r_right,r,r_left + r_right);

            if (diff > 1.0E-8)
            {
                printf("[-] ERROR! Bifurcation rule was not valid !\n");
                return false;
            }
        }
    }
    return true;
}

bool check_null_segments (std::vector<Segment*> s_list)
{
    bool result = false;
    for (uint32_t i = 0; i < s_list.size(); i++)
    {
        if (s_list[i]->length == 0.0)
        {
            printf("[!] WARNING! Segment '%u' has a null length!\n",s_list[i]->id);
            result = true;
        } 
    }
    return result;
}

bool check_collisions_and_fill_feasible_segments (std::vector<Segment*> s_list, Point *p, std::vector<Segment*> &feasible_segments)
{
    std::vector< std::pair<double,uint32_t> > arr;
    for (uint32_t i = 0; i < s_list.size(); i++)
    {
        Segment *s = s_list[i];

        //if (!has_collision(s,p) && s->is_touchable())
        if (!has_collision(s_list,s,p))
        {
            double middle_pos[3];
            s->calc_middle_point(middle_pos);

            double dist = euclidean_norm(middle_pos[0],middle_pos[1],middle_pos[2],\
                                    p->x,p->y,p->z);

            arr.push_back(std::make_pair(dist,s->id));

            feasible_segments.push_back(s);
        }      
    }
    
    // Sort the segments by the distance to the middle point
    std::sort(arr.begin(),arr.end());

    // Get only the closest NCONN segments
    for (uint32_t i = 0; i < arr.size() && feasible_segments.size() < NCONN; i++)
    {
        uint32_t id = arr[i].second;
        
        Segment *s = s_list[id];

        feasible_segments.push_back(s);
    }

    return (feasible_segments.size() == 0) ? false : true;
}

bool has_collision (std::vector<Segment*> s_list, Segment *s, Point *p)
{
    double middle_pos[3];
    s->calc_middle_point(middle_pos);

    for (uint32_t i = 0; i < s_list.size(); i++)
    {
        Segment *cur_segment = s_list[i];
        uint32_t cur_id = cur_segment->id;

        if (s->id != cur_id)
        {
            Point *src = cur_segment->src;
            Point *dest = cur_segment->dest;

            bool intersect = collision_detection(src->x,src->y,src->z,\
                                            dest->x,dest->y,dest->z,\
                                            middle_pos[0],middle_pos[1],middle_pos[2],\
                                            p->x,p->y,p->z);
            if (intersect)
            {
                printf("\t[-] ERROR! Intersection with segment %d !\n",cur_id);
                return true;
            }
        }
    }

    return false;
}

bool connection_search (std::vector<Segment*> s_list, Point *p, const double d_threash)
{
    // Copy the coordinates of the new terminal to an array
    double pos[3];
    pos[0] = p->x;
    pos[1] = p->y;
    pos[2] = p->z;

    for (uint32_t i = 0; i < s_list.size(); i++)
    {
        Segment *s = s_list[i];

        if (!distance_criterion(s,pos,d_threash))
            return false;
    }
    return true;
}

bool distance_criterion (Segment *s, const double pos[], const double d_threash)
{
    double d_proj = calc_dproj(s,pos);

    double d_crit;
    if (d_proj >= 0 && d_proj <= 1)
        d_crit = calc_dortho(s,pos);
    else
        d_crit = calc_dend(s,pos);

    return (d_crit < d_threash) ? false : true;
}

void update_segment_pointers (Segment *iconn, Segment *ibiff, Segment *inew, const bool is_restore)
{
    // Update pointers when building a new segment
    if (!is_restore)
    {
        // iconn:
        iconn->parent = ibiff;
        iconn->src = inew->src;

        // ibiff:
        ibiff->left = inew;     // CONVENTION: Left will always point to terminal
        ibiff->right = iconn;   // CONVENTION: Right will always point to subtree
        if (ibiff->parent != NULL)
        {
            Segment *ibiff_parent = ibiff->parent;
            if (ibiff_parent->left->id == iconn->id)
                ibiff_parent->left = ibiff;
            if (ibiff_parent->right->id == iconn->id)
                ibiff_parent->right = ibiff;
        }
    }
    // Update pointers when restoring the tree state
    else
    {
        Segment *ibiff_par = ibiff->parent;

        // iconn:
        iconn->parent = ibiff->parent;
        iconn->src = ibiff->src;
        iconn->beta = ibiff->beta;

        // ibiff_par:
        if (ibiff_par != NULL)
        {
            if (ibiff_par->right == ibiff)
                ibiff_par->right = iconn;
            if (ibiff_par->left == ibiff)
                ibiff_par->left = iconn;
        }
    }
}

void get_segment_length (std::vector<Segment*> s_list, std::vector<double> &segments)
{
    for (uint32_t i = 0; i < s_list.size(); i++)
    {
        Segment *cur_segment = s_list[i];

        // Segment length (This value is given in {m})
        double length = cur_segment->length;
        if (length == 0.0)
        {
            printf("[!] WARNING! Segment length is equal zero! Indexes: %u %u\n",cur_segment->src->id,cur_segment->dest->id);
        }

        segments.push_back(length * M_TO_MM);
    }
}

void get_bifurcation_angles(std::vector<Segment*> s_list, std::vector<double> &angles)
{
    for (uint32_t i = 0; i < s_list.size(); i++)
    {
        Segment *cur_segment = s_list[i];

        if (cur_segment->left != NULL && cur_segment->right != NULL)
        {
            double u[3], v[3], angle;
            cur_segment->left->calc_unitary_vector(u);
            cur_segment->right->calc_unitary_vector(v);
            angle = calc_angle_between_vectors(u,v);

            angles.push_back(angle);
        }
    }
}

Point* generate_bifurcation_node (std::vector<Point*> p_list, Segment *iconn, LocalOptimizationConfig *local_opt_config, const bool using_local_optimization)
{
    uint32_t cur_num_nodes = p_list.size(); 
    Point *result = new Point(cur_num_nodes);

    if (using_local_optimization)
    {
        // The local optimization was not executed yet
        if (local_opt_config->first_call)
        {
            double middle_pos[3];
            iconn->calc_middle_point(middle_pos);
            result->setCoordinate(middle_pos);
        }
        // The best bifurcation position is already stored in the 'best_pos' array
        else
        {
            double *best_pos = local_opt_config->best_pos;
            result->setCoordinate(best_pos);

            // Reset the 'first_call' flag for the next point
            local_opt_config->first_call = true;
        }
    }
    else
    {
        double middle_pos[3];
        iconn->calc_middle_point(middle_pos);
        result->setCoordinate(middle_pos);
    }
    result->setLAT(iconn->calc_middle_point_lat());
    result->setActive(false);

    return result;
}

Point* generate_terminal_node (std::vector<Point*> p_list, Point *p)
{
    uint32_t cur_num_nodes = p_list.size();
    Point *result = new Point(cur_num_nodes);

    double pos[3];
    pos[0] = p->x;
    pos[1] = p->y;
    pos[2] = p->z;

    result->setCoordinate(pos);
    result->setActive(p->is_active);
    result->setLAT(p->lat);

    return result;
}

Segment* get_terminal (std::vector<Point*> p_list, std::vector<Segment*> s_list, const double pos[])
{
    uint32_t min_id = 0;
    double min_dist = __DBL_MAX__;

    // Find the index of the closest point to the one given
    for (uint32_t i = 0; i < p_list.size(); i++)
    {
        double middle_pos[3];
        middle_pos[0] = p_list[i]->x*M_TO_UM;
        middle_pos[1] = p_list[i]->y*M_TO_UM;
        middle_pos[2] = p_list[i]->z*M_TO_UM;

        double dist = euclidean_norm(middle_pos[0],middle_pos[1],middle_pos[2],\
                                    pos[0],pos[1],pos[2]);
        
        if (dist < min_dist)
        {
            min_dist = dist;
            min_id = i;
        }
    }

    for (uint32_t i = 0; i < s_list.size(); i++)
    {
        if (s_list[i]->dest->id == min_id)
            return s_list[i];
    }
    for (uint32_t i = 0; i < s_list.size(); i++)
    {
        if (s_list[i]->src->id == min_id)
            return s_list[i];
    }
    return NULL;
}

uint32_t update_ndist (Segment *s, const bool is_restore)
{
    // Update pointers when building a new segment
    if (!is_restore)
    {
        Segment *tmp = s;
        while (tmp->parent != NULL)
        {
            tmp->ndist++;
            tmp = tmp->parent;
        }
        tmp->ndist++;

        // Update the current number of terminals in the network
        return tmp->ndist;
    }
    // Update pointers when restoring the tree state
    else
    {
        Segment *tmp = s->parent;
        if (tmp != NULL)
        {
            while (tmp->parent != NULL)
            {
                tmp->ndist--;
                tmp = tmp->parent;
            }
            tmp->ndist--;
            return tmp->ndist;
        }
        else
        {
            return s->ndist;
        }
    }
}