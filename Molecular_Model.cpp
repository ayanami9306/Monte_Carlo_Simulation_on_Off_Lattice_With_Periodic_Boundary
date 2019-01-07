#include "Molecular_Model.hpp"
#include <sys/time.h>
#include <unistd.h>

Model_Segment::Model_Segment()
{
    inv_sin_MAX = 1.0 / 65536.0;
    exp_divide = 1.0 / 1024.0;
    Half_Boundary_X = BOUNDARY_X/2.0;
    Half_Boundary_Y = BOUNDARY_Y/2.0;
    Half_Boundary_Z = BOUNDARY_Z/2.0;
    srand((unsigned int)time(NULL));
}

Model_Segment::~Model_Segment()
{
    
}

