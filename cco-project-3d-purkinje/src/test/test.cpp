#include "test.h"

#include "../cco/cco.h"

void test1 (CCO_Network *the_network)
{
    double pos1[3] = {0,0,0};
    double pos2[3] = {-3,-3,0};
    double pos3[3] = {-1,-1,0};
    double pos4[3] = {1,-3,0};
    double pos5[3] = {-2,-2,0};
    double pos6[3] = {-2,-4,0};
    double pos7[3] = {0,-2,0};
    double pos8[3] = {-1,-3,0};

    Point *A = new Point(0,pos1);
    Point *B = new Point(1,pos2);
    Point *C = new Point(2,pos3);
    Point *D = new Point(3,pos4);
    Point *E = new Point(4,pos5);
    Point *F = new Point(5,pos6);
    Point *G = new Point(6,pos7);
    Point *H = new Point(7,pos8);

    the_network->point_list.push_back(A);
    the_network->point_list.push_back(B);
    the_network->point_list.push_back(C);
    the_network->point_list.push_back(D);
    the_network->point_list.push_back(E);
    the_network->point_list.push_back(F);
    the_network->point_list.push_back(G);
    the_network->point_list.push_back(H);

    Segment *s1 = new Segment(0,A,C,NULL,NULL,NULL);
    Segment *s2 = new Segment(1,C,E,NULL,NULL,NULL);
    Segment *s3 = new Segment(2,C,G,NULL,NULL,NULL);
    Segment *s4 = new Segment(3,G,H,NULL,NULL,NULL);
    Segment *s5 = new Segment(4,G,D,NULL,NULL,NULL);
    Segment *s6 = new Segment(5,E,B,NULL,NULL,NULL);
    Segment *s7 = new Segment(6,E,F,NULL,NULL,NULL);

    s1->parent = NULL;
    s1->left = s2;
    s1->right = s3;

    s2->parent = s1;
    s2->left = s6;
    s2->right = s7;

    s3->parent = s1;
    s3->left = s4;
    s3->right = s5;

    s4->parent = s3;
    s4->left = NULL;
    s4->right = NULL;

    s5->parent = s3;
    s5->left = NULL;
    s5->right = NULL;

    s6->parent = s2;
    s6->left = NULL;
    s6->right = NULL;

    s7->parent = s2;
    s7->left = NULL;
    s7->right = NULL;

    the_network->segment_list.push_back(s1);
    the_network->segment_list.push_back(s2);
    the_network->segment_list.push_back(s3);
    the_network->segment_list.push_back(s4);
    the_network->segment_list.push_back(s5);
    the_network->segment_list.push_back(s6);
    the_network->segment_list.push_back(s7);

    the_network->print();
}