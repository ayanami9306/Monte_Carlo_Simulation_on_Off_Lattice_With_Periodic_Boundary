#include "Molecular_Model.hpp"

Model_Segment System;

int main(int argc, char*argv[])
{
    char Input_Param_File[50];
    strcpy(Input_Param_File, argv[1]);
    System.Input_Params(Input_Param_File);
    System.Initialize_System(1);
    System.MonteCarlo();
}
