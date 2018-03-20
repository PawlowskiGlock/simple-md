#ifndef SYSTEM_HPP
#define SYSTEM_HPP

// INCLUDE ---------------------------------------------------------------------
#include <vector>


// NAMESPACE -------------------------------------------------------------------
namespace exercises {
    // CONSTANTS
    const double M = 15.035; // mass of CH3 fragment

    const double K = 3e5; // bond constant and length
    const double B0 = 0.153;

    const double C12 = 33.570436e-6; // LJ-constants
    const double C6 = 0.0099161764;

    const int DIM = 3; // dimensionality

    /* INFO
     * The classes Atom, Bond and Molecule are not full-blown classes but
     * rather structs with a constructor. Their attributes need to be modified
     * by the System class, which is therefore declared as a friend.
     */

    // ATOM --------------------------------------------------------------------
    class Atom {
    private:
        double mass;
        std::vector<double> pos;
        std::vector<double> vel;
        std::vector<double> force;
    public:
        Atom(double mass_, std::vector<double> pos_, std::vector<double> vel_);
        friend class System;
    };


    // BOND --------------------------------------------------------------------
    class Bond {
    private:
        double k;
        double b0;
        Atom *atom1;
        Atom *atom2;
    public:
        Bond(double k_, double b0_, Atom *atom1_, Atom *atom2_);
        friend class System;
    };


    // MOLECULE ----------------------------------------------------------------
    class Molecule {
    private:
        std::vector<Atom*> atoms;
        std::vector<Bond*> bonds;
    public:
        Molecule(std::vector<Atom*> atoms_, std::vector<Bond*> bonds_);
        friend class System;
    };


    // SYSTEM ------------------------------------------------------------------
    class System {
    private:
        std::vector<Molecule> mols;
        void bond_force(Bond *bond);
        void lj_force(Atom *atom1, Atom *atom2);
        double bond_energy_x2(Bond *bond);
        double lj_energy(Atom *atom1, Atom *atom2);
        void kick(double dt);
        void drift(double dt);
        double kinetic_energy();
        double potential_energy();
        void print_atoms(int step, double timestep);
    public:
        System(std::vector<Molecule> mols_);
        void md(int nsteps, double timestep);
    };
}

#endif // SYSTEM_HPP
