
#include "Molecular_Model.hpp"

static double sin90_table[] = {
    0,  1143,  2287,  3429,  4571,  5711,  6850,  7986,  9120, 10252,
    11380, 12504, 13625, 14742, 15854, 16961, 18064, 19160, 20251, 21336,
    22414, 23486, 24550, 25606, 26655, 27696, 28729, 29752, 30767, 31772,
    32767, 33753, 34728, 35693, 36647, 37589, 38521, 39440, 40347, 41243,
    42125, 42995, 43852, 44695, 45525, 46340, 47142, 47929, 48702, 49460,
    50203, 50931, 51643, 52339, 53019, 53683, 54331, 54963, 55577, 56175,
    56755, 57319, 57864, 58393, 58903, 59395, 59870, 60326, 60763, 61183,
    61583, 61965, 62328, 62672, 62997, 63302, 63589, 63856, 64103, 64331,
    64540, 64729, 64898, 65047, 65176, 65286, 65376, 65446, 65496, 65526,
    65536
};

static int neighbor_cell_offset[][3] =
{{0, 0, 0}, {0, 0, 1}, {0, 0, -1}, {0, 1, 0}, {0, -1, 0}, {1, 0, 0}, {-1, 0, 0},
    {0, 1, 1}, {0, -1, 1}, {0, 1, -1}, {0, -1, -1},
    {1, 0, 1}, {-1, 0, 1}, {1, 0, -1}, {-1, 0, -1},
    {1, 1, 0}, {1, -1, 0}, {-1, 1, 0}, {-1, -1, 0},
    {1, 1, 1}, {1, 1, -1}, {1, -1, 1}, {-1, 1, 1},
    {-1, -1, 1}, {-1, 1,-1}, {1, -1, -1}, {-1, -1, -1}
};

double Model_Segment::optimize_cos(int angle)
{
    if(angle <= 90) return sin90_table[90-angle] * inv_sin_MAX;
    if(angle <= 180) return -sin90_table[angle-90] * inv_sin_MAX;
    if(angle <= 270) return -sin90_table[270-angle] * inv_sin_MAX;
    return sin90_table[angle-270] * inv_sin_MAX;
}

double Model_Segment::optimize_sin(int angle)
{
    if(angle <= 90) return sin90_table[angle] * inv_sin_MAX;
    if(angle <= 180) return sin90_table[180-angle] * inv_sin_MAX;
    if(angle <= 270) return -sin90_table[angle-180] * inv_sin_MAX;
    return -sin90_table[360-angle] * inv_sin_MAX;
}

double Model_Segment::fast_exp(double x)
{
    x = 1.0 + x * exp_divide;
    x *= x;
    x *= x;
    x *= x;
    x *= x;
    x *= x;
    x *= x;
    x *= x;
    x *= x;
    x *= x;
    x *= x;
    return x;
}

void Model_Segment::MC_CONSTRUCT_CELL()
{
    cellnum[0] = (int)(BOUNDARY_X * inv_rcut);
    cellnum[1] = (int)(BOUNDARY_Y * inv_rcut);
    cellnum[2] = (int)(BOUNDARY_Z * inv_rcut);
    
    int total_cell_num = cellnum[0] * cellnum[1] * cellnum[2];
    for(int i=0; i<total_cell_num; i++)
    {
        list<int> * temp = new list<int>;
        celllist.push_back(temp);
    }

    int cellposition[3];
    for(int i=0; i<nParticle;i++)
    {
        //X, Y, Z
        cellposition[0] = (int)(Segment[i].coordinate[0] * inv_rcut);
        cellposition[1] = (int)(Segment[i].coordinate[1] * inv_rcut);
        cellposition[2] = (int)(Segment[i].coordinate[2] * inv_rcut);
        celllist[cellposition[2] * cellnum[0] * cellnum[1] + cellposition[1] * cellnum[0] + cellposition[0]]->push_back(i);
    }
}

double Model_Segment::MC_NEIGHBORPOTENTIAL(int index)
{
    int cellposition[3]{
        (int)(Segment[index].coordinate[0] * inv_rcut), (int)(Segment[index].coordinate[1] * inv_rcut), (int)(Segment[index].coordinate[2] * inv_rcut)
    };
    double distance2, pot = 0;
    int temp_cellposition[3], temp_cellindex;
    list<int> * temp_cell;
    list<int>::iterator itor, itor_end;
    
    //i = 0(self cell)
    temp_cellposition[0] = cellposition[0] + neighbor_cell_offset[0][0];
    temp_cellposition[1] = cellposition[1] + neighbor_cell_offset[0][1];
    temp_cellposition[2] = cellposition[2] + neighbor_cell_offset[0][2];
    
    if(temp_cellposition[0] >= cellnum[0]) temp_cellposition[0] -= cellnum[0];
    if(temp_cellposition[0] < 0) temp_cellposition[0] += cellnum[0];
    if(temp_cellposition[1] >= cellnum[1]) temp_cellposition[1] -= cellnum[1];
    if(temp_cellposition[1] < 0) temp_cellposition[1] += cellnum[1];
    if(temp_cellposition[2] >= cellnum[2]) temp_cellposition[2] -= cellnum[2];
    if(temp_cellposition[2] < 0) temp_cellposition[2] += cellnum[2];
    
    temp_cellindex = temp_cellposition[2] * cellnum[0] * cellnum[1] + temp_cellposition[1] * cellnum[0] + temp_cellposition[0];
    temp_cell = celllist[temp_cellindex];
    itor_end = temp_cell->end();
    for(itor = temp_cell->begin(); itor != itor_end; itor++)
    {
        if(*itor != index)
        {
            distance2 = Get_Distance2(index, *itor);
            if(distance2 <= rcut2)
            {
                double inv_distance2 = 1.0 / distance2;
                double inv_distance6 = pow(inv_distance2, 3.0);
                pot += 4 * epsilon * inv_distance6 * (inv_distance6 - 1.0) - potential_rcut;
            }
        }
    }
    
    //i != 0
    for(int i=1; i<27; i++)
    {
        temp_cellposition[0] = cellposition[0] + neighbor_cell_offset[i][0];
        temp_cellposition[1] = cellposition[1] + neighbor_cell_offset[i][1];
        temp_cellposition[2] = cellposition[2] + neighbor_cell_offset[i][2];
        
        if(temp_cellposition[0] >= cellnum[0]) temp_cellposition[0] -= cellnum[0];
        if(temp_cellposition[0] < 0) temp_cellposition[0] += cellnum[0];
        if(temp_cellposition[1] >= cellnum[1]) temp_cellposition[1] -= cellnum[1];
        if(temp_cellposition[1] < 0) temp_cellposition[1] += cellnum[1];
        if(temp_cellposition[2] >= cellnum[2]) temp_cellposition[2] -= cellnum[2];
        if(temp_cellposition[2] < 0) temp_cellposition[2] += cellnum[2];

        
        temp_cellindex = temp_cellposition[2] * cellnum[0] * cellnum[1] + temp_cellposition[1] * cellnum[0] + temp_cellposition[0];
        
        temp_cell = celllist[temp_cellindex];
        itor_end = temp_cell->end();
        for(itor = temp_cell->begin(); itor != itor_end; itor++)
        {
            distance2 = Get_Distance2(index, *itor);
            if(distance2 <= rcut2)
            {
                double inv_distance2 = 1.0 / distance2;
                double inv_distance6 = pow(inv_distance2, 3.0);
                pot += 4 * epsilon * inv_distance6 * (inv_distance6 - 1.0) - potential_rcut;
            }
        }
    }
    
    return pot;
}
void Model_Segment::MonteCarlo()
{
    mt19937 engine((unsigned int)time(NULL));
    uniform_int_distribution<int> distribution(0, nParticle - 1);
    auto generator = bind(distribution, engine);
    pot_step = 0;
    Set_Params();
    //Mol2_File_Write(true);
    double * temp_coordinate = (double *)malloc(3 * sizeof(double));
    int num_prev, rSegment;
    double distance_prev2;
    double pot_pre, pot_after, harmonic_pre, harmonic_after, theta, phi, radius, sin_theta;
    int selected_linked_segment_num;
    int pre_cell_position[3], after_cell_position[3];
    Node *temp_Segment;
    list<int> * pre_cell, * after_cell;
    list<int>::iterator itor;
    int status = 0;
    MC_CONSTRUCT_CELL();
    for(int mccount=0; mccount<=Limit_Cycle; mccount++)
    {
        for(int i=0; i<nParticle; i++)
        {
            //calc before potential
            pot_pre = 0; pot_after = 0; harmonic_pre = 0; harmonic_after = 0;
            rSegment = generator();
            temp_Segment = &Segment[rSegment];
            //lennar -jones potential
            pot_pre += MC_NEIGHBORPOTENTIAL(rSegment);
            
            //harmonic potential
            selected_linked_segment_num = temp_Segment->linked_segment_num;
            for(int j=0; j<selected_linked_segment_num ;j++)
            {
                num_prev = temp_Segment->linked_segment[j];
                //ba vector
                distance_prev2 = Get_Distance2(num_prev, rSegment);
                harmonic_pre += distance_prev2 - 2 * sqrt(distance_prev2 * bond_length_harmonic2);
            }
            pot_pre += k_harmonic * (harmonic_pre + selected_linked_segment_num * bond_length_harmonic2);
            
            //random hopping
            theta = (double)generator() * inv_RAND_MAX_scaled_by_360;
            phi = (double)generator() * inv_RAND_MAX_scaled_by_360;
            radius = (double)generator() * inv_RAND_MAX_scaled_by_rMax;
            memcpy(temp_coordinate, temp_Segment->coordinate, 3 * sizeof(double));
            sin_theta = optimize_sin(theta);
            temp_Segment->coordinate[0] += radius * sin_theta * optimize_cos(phi);
            temp_Segment->coordinate[1] += radius * sin_theta * optimize_sin(phi);
            temp_Segment->coordinate[2] += radius * optimize_cos(theta);
            
            //Apply periodic boundary
            if(temp_Segment->coordinate[0] >= BOUNDARY_X) temp_Segment->coordinate[0] -= BOUNDARY_X;
            if(temp_Segment->coordinate[0] < 0) temp_Segment->coordinate[0] += BOUNDARY_X;
            if(temp_Segment->coordinate[1] >= BOUNDARY_Y) temp_Segment->coordinate[1] -= BOUNDARY_Y;
            if(temp_Segment->coordinate[1] < 0) temp_Segment->coordinate[1] += BOUNDARY_Y;
            if(temp_Segment->coordinate[2] >= BOUNDARY_Z) temp_Segment->coordinate[2] -= BOUNDARY_Z;
            if(temp_Segment->coordinate[2] < 0) temp_Segment->coordinate[2] += BOUNDARY_Z;
            
            //cell update
            pre_cell_position[0] = (int)(temp_coordinate[0] * inv_rcut);
            pre_cell_position[1] = (int)(temp_coordinate[1] * inv_rcut);
            pre_cell_position[2] = (int)(temp_coordinate[2] * inv_rcut);
            after_cell_position[0] = (int)(temp_Segment->coordinate[0] * inv_rcut);
            after_cell_position[1] = (int)(temp_Segment->coordinate[1] * inv_rcut);
            after_cell_position[2] = (int)(temp_Segment->coordinate[2] * inv_rcut);
            
            if(pre_cell_position[0] != after_cell_position[0] || pre_cell_position[1] != after_cell_position[1] || pre_cell_position[2] != after_cell_position[2])
            {
                pre_cell = celllist[pre_cell_position[2] * cellnum[0] * cellnum[1] + pre_cell_position[1] * cellnum[0] + pre_cell_position[0]];
                after_cell = celllist[after_cell_position[2] * cellnum[0] * cellnum[1] + after_cell_position[1] * cellnum[0] + after_cell_position[0]];
                
                itor = find(pre_cell->begin(), pre_cell->end(), rSegment);
                
                pre_cell->erase(itor);
                after_cell->push_back(rSegment);
                
                status = 1;
            }
            
            //calc after potential
            pot_after += MC_NEIGHBORPOTENTIAL(rSegment);
            for(int j=0; j<selected_linked_segment_num; j++)
            {
                int num_prev = temp_Segment->linked_segment[j];
                //ba vector
                double distance_prev2 = Get_Distance2(num_prev, rSegment);
                harmonic_after += distance_prev2 - 2 * sqrt(distance_prev2 * bond_length_harmonic2);
            }
            pot_after += k_harmonic * (harmonic_after + selected_linked_segment_num * bond_length_harmonic2);
            
            double delta_pot = (-pot_after + pot_pre) * inv_kT_0;
            
            //reject
            if(delta_pot <= 0)
            {
                if(exp(delta_pot) < (double)generator()*inv_RAND_MAX)
                {
                    memcpy(temp_Segment->coordinate, temp_coordinate, 3 * sizeof(double));
                    pot_after = pot_pre;
                    
                    if(status)
                    {
                        pre_cell->push_back(rSegment);
                        after_cell->pop_back();
                    }
                }
            }
            pot_step += pot_after;
            status = 0;
            
        }
        if(!(mccount % step_AVG) && mccount)
        {
            Write_State(mccount);
            //Mol2_File_Write(false);
            Evaluate_Properties(mccount);
        }
    }
    free(temp_coordinate);
}
