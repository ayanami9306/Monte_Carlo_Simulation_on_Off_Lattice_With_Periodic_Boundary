/***********************************************
 
 LENGTH IS PROPORTIONAL TO SIGMA
 ENERGY IS PROPORTIONAL TO EPSILON
 
 **********************************************/
#ifndef Molecular_Model_hpp
#define Molecular_Model_hpp

#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <list>
#include <vector>
#include <algorithm>
#include <random>

#define NUM_MOLECULE 100000
#define PI 3.141592
#define double_PI 6.283184
#define BOUNDARY_X 100
#define BOUNDARY_Y 100
#define BOUNDARY_Z 100

using namespace std;

#endif /* Molecular_Model_hpp */

typedef struct Node
{
    double coordinate[3];
    int linked_segment_num;
    int linked_segment[10];
    int segment_type;
}Node;

class Model_Segment
{
    //variables
private:
    Node Segment[NUM_MOLECULE];
    int nParticle, step_AVG, Limit_Cycle, dp_backbone, dp_dendron, generation_dendron, number_branch, space_sidechain, addi_frag;
    double rcut, epsilon, rMax, inv_RAND_MAX, inv_RAND_MAX_scaled_by_360, inv_RAND_MAX_scaled_by_rMax, inv_rcut;
    double rcut2, inv_step_AVG, potential_rcut, kT_0, inv_kT_0, pot_step, time_Now, bond_length, bond_length_harmonic, k_harmonic, bond_length_harmonic2, inv_sin_MAX, exp_divide, Half_Boundary_X, Half_Boundary_Y, Half_Boundary_Z;
    char write_filename[100], dat_filename[100], traj_filename[100];
    int cellnum[3];
    vector<list<int> *> celllist;
    
    //methods
public:
    Model_Segment();
    ~Model_Segment();
    void Initialize_System(int nTypes);
    void Input_Params(char *filename);
    void Read_State(char *filename);
    void MonteCarlo();
    
private:
    void tolower(char *data);
    void Set_Params();
    void Write_State(int Cycle);
    double Get_Distance2(int num1, int num2);
    void recursive_branch(int parent, int generation);
    void Evaluate_Properties(int mccount);
    void Mol2_File_Write(bool is_new_file);
    void Mol2_File_Read(char *filename);
    double MC_NEIGHBORPOTENTIAL(int index);
    void Periodic_Boundary(double *coordinate);
    void MC_CONSTRUCT_CELL();
    double optimize_sin(int angle);
    double optimize_cos(int angle);
    double fast_exp(double x);
};
