#include "Molecular_Model.hpp"

void Model_Segment::Input_Params(char *filename)
{
    FILE *fp = fopen(filename, "r");
    char property_value[50];
    char property_name[50];
    while(fscanf(fp, "%s %s\n", property_name, property_value) > 0)
    {
        tolower(property_name);
        if(!strcmp(property_name, "rcut")) rcut = atof(property_value);
        else if(!strcmp(property_name, "step_avg")) step_AVG = atoi(property_value);
        else if(!strcmp(property_name, "limit_cycle")) Limit_Cycle = atof(property_value);
        else if(!strcmp(property_name, "write_filename")) strcpy(write_filename, property_value);
        else if(!strcmp(property_name, "bond_length_harmonic")) bond_length_harmonic = atof(property_value);
        else if(!strcmp(property_name, "k_harmonic")) k_harmonic = atof(property_value);
        else if(!strcmp(property_name, "kt_0")) kT_0 = atof(property_value);
        else if(!strcmp(property_name, "dp_backbone")) dp_backbone = atoi(property_value);
        else if(!strcmp(property_name, "dp_dendron")) dp_dendron = atoi(property_value);
        else if(!strcmp(property_name, "generation_dendron")) generation_dendron = atoi(property_value);
        else if(!strcmp(property_name, "number_branch")) number_branch   = atoi(property_value);
        else if(!strcmp(property_name, "space_sidechain")) space_sidechain   = atoi(property_value);
        else if(!strcmp(property_name, "epsilon")) epsilon = atof(property_value);
        else if(!strcmp(property_name, "addi_frag")) addi_frag = atoi(property_value);
    }
    fclose(fp);
}

void Model_Segment::Set_Params()
{
    rcut2 = rcut*rcut; //rcut unit : sigma
    rMax = bond_length_harmonic;
    inv_rcut = 1.0 / rcut;
    double inverse_rcut6 = 1.0 / pow(rcut2, 3.0);
    potential_rcut = 4 * epsilon * inverse_rcut6 * (inverse_rcut6 - 1.0);
    
    bond_length_harmonic2 = bond_length_harmonic * bond_length_harmonic;
    inv_RAND_MAX_scaled_by_360 = 360.0 / (double)RANDOM_MAX;
    inv_RAND_MAX_scaled_by_rMax = rMax / (double)RANDOM_MAX;
    inv_RAND_MAX = 1.0 / (double)RANDOM_MAX;
    
    inv_step_AVG = 1.0 / (double)step_AVG;
    inv_kT_0 = 1.0 / kT_0;
    
    sprintf(dat_filename, "%s_DAT.csv", write_filename);
    sprintf(traj_filename, "%s.mol2", write_filename);
    if(time_Now == 0)
    {
        FILE *fp = fopen(dat_filename, "w");
        fprintf(fp, "TIME,POT\n");
        fclose(fp);
    }
}
