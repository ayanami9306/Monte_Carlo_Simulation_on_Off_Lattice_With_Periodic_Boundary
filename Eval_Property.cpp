#include "Molecular_Model.hpp"

void Model_Segment::Evaluate_Properties(int mccount)
{
    pot_step = (pot_step * inv_step_AVG)/nParticle;
    FILE *fp = fopen(dat_filename, "a");
    fprintf(fp, "%d,%.4lf\n", mccount,  pot_step);
    fclose(fp);
    
    pot_step = 0;
}

double Model_Segment::Get_Distance2(int num1, int num2)
{
    double dr[3];
    dr[0] = fabs(Segment[num1].coordinate[0] - Segment[num2].coordinate[0]);
    dr[1] = fabs(Segment[num1].coordinate[1] - Segment[num2].coordinate[1]);
    dr[2] = fabs(Segment[num1].coordinate[2] - Segment[num2].coordinate[2]);
    if(dr[0] > Half_Boundary_X) dr[0] = BOUNDARY_X - dr[0];
    if(dr[1] > Half_Boundary_Y) dr[1] = BOUNDARY_Y - dr[1];
    if(dr[2] > Half_Boundary_Z) dr[2] = BOUNDARY_Z - dr[2];
    return dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];
}
