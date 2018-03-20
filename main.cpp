// INCLUDE ---------------------------------------------------------------------
#include <iostream>
#include "system.hpp"


// USING -----------------------------------------------------------------------
using namespace std;
using namespace exercises;


// MAIN ------------------------------------------------------------------------
int main()
{
    // scenario 1
    vector<double> pos1 = { 0.0000,  0.0000,  0.0000};
    vector<double> pos2 = { 0.1530,  0.0000,  0.0000};
    vector<double> pos3 = { 2.0000,  0.0000,  0.0000};
    vector<double> pos4 = { 2.1530,  0.0000,  0.0000};

    vector<double> vel1 = { 0.0000,  0.0000,  0.0000};
    vector<double> vel2 = { 0.0000,  0.0000,  0.0000};
    vector<double> vel3 = {-1.0000,  0.0000,  0.0000};
    vector<double> vel4 = {-1.0000,  0.0000,  0.0000};

    /*
     // scenario 2
     vector<double> pos1 = { 0.0000, -0.0765,  0.0000};
     vector<double> pos2 = { 0.0000,  0.0765,  0.0000};
     vector<double> pos3 = { 2.0000,  0.0000,  0.0000};
     vector<double> pos4 = { 2.1530,  0.0000,  0.0000};

     vector<double> vel1 = { 0.0000,  0.0000,  0.0000};
     vector<double> vel2 = { 0.0000,  0.0000,  0.0000};
     vector<double> vel3 = {-1.0000,  0.0000,  0.0000};
     vector<double> vel4 = {-1.0000,  0.0000,  0.0000};

     // scenario 3
     vector<double> pos1 = { 0.0000,  0.0000,  0.0000};
     vector<double> pos2 = { 0.0000,  0.1530,  0.0000};
     vector<double> pos3 = { 2.0000,  0.0000,  0.0000};
     vector<double> pos4 = { 2.1530,  0.0000,  0.0000};

     vector<double> vel1 = { 0.0000,  0.0000,  0.0000};
     vector<double> vel2 = { 0.0000,  0.0000,  0.0000};
     vector<double> vel3 = {-1.0000,  0.0000,  0.0000};
     vector<double> vel4 = {-1.0000,  0.0000,  0.0000};
     */

    // create Atoms
    Atom a1(M, pos1, vel1);
    Atom a2(M, pos2, vel2);
    Atom a3(M, pos3, vel3);
    Atom a4(M, pos4, vel4);

    // create Bonds
    Bond b12(K, B0, &a1, &a2);
    Bond b34(K, B0, &a3, &a4);

    // create Molecules
    vector<Atom*> atoms_m1 = {&a1, &a2};
    vector<Atom*> atoms_m2 = {&a3, &a4};
    vector<Bond*> bonds_m1 = {&b12};
    vector<Bond*> bonds_m2 = {&b34};

    Molecule m1(atoms_m1, bonds_m1);
    Molecule m2(atoms_m2, bonds_m2);

    // create System
    vector<Molecule> mols_of_sys = {m1, m2};

    System sys(mols_of_sys);
    /* Now all the vectors of Atom* and Bond* as well as the Molecules could
     * be deleted, the System has copied the molecules which hold pointers
     * to the original Atom and Bond objects
     */
    
    // run simulation
    sys.md(500, 0.002);
    
    return 0;
}
