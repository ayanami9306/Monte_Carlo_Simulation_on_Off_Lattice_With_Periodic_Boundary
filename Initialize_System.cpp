#include "Molecular_Model.hpp"

void Model_Segment::Initialize_System(int nTypes)
{
    nParticle = 0;
    inv_RAND_MAX = 1.0 / (double)RAND_MAX;
    //dedronized polymer
    if(nTypes == 1)
    {
        double x = Half_Boundary_X, y = Half_Boundary_Y, z = Half_Boundary_Z;
        for(int i=0;i<dp_backbone;i++)
        {
            double torsional_angle = ((double)rand()*inv_RAND_MAX)*2*PI;
            double bond_angle = ((double)rand()*inv_RAND_MAX)*2*PI;
            Segment[i].coordinate[0] = x + bond_length_harmonic * sin(bond_angle)*cos(torsional_angle);
            Segment[i].coordinate[1] = y + bond_length_harmonic * sin(bond_angle)*sin(torsional_angle);
            Segment[i].coordinate[2] = z + bond_length_harmonic * cos(bond_angle);
            Segment[i].linked_segment_num = 0;
            Segment[i].segment_type = 0;
            
            //Apply periodic boundary
            if(Segment[i].coordinate[0] >= BOUNDARY_X) Segment[i].coordinate[0] -= BOUNDARY_X;
            if(Segment[i].coordinate[0] < 0) Segment[i].coordinate[0] += BOUNDARY_X;
            if(Segment[i].coordinate[1] >= BOUNDARY_Y) Segment[i].coordinate[1] -= BOUNDARY_X;
            if(Segment[i].coordinate[1] < 0) Segment[i].coordinate[1] += BOUNDARY_X;
            if(Segment[i].coordinate[2] >= BOUNDARY_Z) Segment[i].coordinate[2] -= BOUNDARY_X;
            if(Segment[i].coordinate[2] < 0) Segment[i].coordinate[2] += BOUNDARY_X;
            
            x = Segment[i].coordinate[0];
            y = Segment[i].coordinate[1];
            z = Segment[i].coordinate[2];
            if(i != dp_backbone-1) Segment[i].linked_segment[Segment[i].linked_segment_num++] = i+1;
            if(i != 0) Segment[i].linked_segment[Segment[i].linked_segment_num++] = i-1;
            nParticle++;
        }
        for(int i=0; i<dp_backbone;i+=space_sidechain) recursive_branch(i, generation_dendron);
    }
    printf("DP : %d\n", (nParticle / dp_backbone) - 1);
}

void Model_Segment::recursive_branch(int parent, int generation)
{
    if(generation)
    {
        int previous_node = parent;
        int dp_node = dp_dendron;
        if(generation == generation_dendron) dp_node = dp_dendron + addi_frag;
        for(int i=0; i<dp_node; i++)
        {
            double torsional_angle = ((double)rand()*inv_RAND_MAX)*2*PI;
            double bond_angle = ((double)rand()*inv_RAND_MAX)*2*PI;
            Segment[nParticle].coordinate[0] = Segment[previous_node].coordinate[0] + bond_length_harmonic * sin(bond_angle)*cos(torsional_angle);
            Segment[nParticle].coordinate[1] = Segment[previous_node].coordinate[1] + bond_length_harmonic * sin(bond_angle)*sin(torsional_angle);
            Segment[nParticle].coordinate[2] = Segment[previous_node].coordinate[2] + bond_length_harmonic * cos(bond_angle);
            Segment[nParticle].linked_segment_num = 0;
            Segment[nParticle].segment_type = 1;
            Segment[previous_node].linked_segment[Segment[previous_node].linked_segment_num++] = nParticle;
            Segment[nParticle].linked_segment[Segment[nParticle].linked_segment_num++] = previous_node;
            previous_node = nParticle;
            nParticle++;
        }
        for(int i=0; i<number_branch; i++) recursive_branch(previous_node, generation-1);
    }
    else return;
}
