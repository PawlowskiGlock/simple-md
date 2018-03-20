// INCLUDE ---------------------------------------------------------------------
#include <iostream>
#include <iomanip>
#include <cmath>

#include "system.hpp"


// USING -----------------------------------------------------------------------
using namespace std;
using namespace exercises;


// ATOM ------------------------------------------------------------------------
exercises::Atom::Atom(double mass_, vector<double> pos_, vector<double> vel_)
: mass(mass_), pos(pos_), vel(vel_), force(3, 0.0) {}


// BOND ------------------------------------------------------------------------
exercises::Bond::Bond(double k_, double b0_, Atom *atom1_, Atom *atom2_)
: k(k_), b0(b0_), atom1(atom1_), atom2(atom2_) {}


// MOLECULE --------------------------------------------------------------------
exercises::Molecule::Molecule(vector<Atom*> atoms_, vector<Bond*> bonds_)
: atoms(atoms_), bonds(bonds_) {}


// SYSTEM ----------------------------------------------------------------------
exercises::System::System(std::vector<Molecule> mols_)
: mols(mols_) {}

void exercises::System::bond_force(Bond *bond) {
    vector<double> r(3);
    double b = 0.0;
    double F = 0.0;
    for (int i = 0; i < DIM; ++i) {
        r[i] = bond->atom1->pos[i] - bond->atom2->pos[i];
        b += r[i] * r[i];
    }
    b = sqrt(b);
    F = bond->k * (b - bond->b0);
    for (int i = 0; i < DIM; ++i) {
        r[i] /= b;
        bond->atom1->force[i] -= F * r[i];
        bond->atom2->force[i] += F * r[i];
    }
    return;
}

void exercises::System::lj_force(Atom *atom1, Atom *atom2) {
    vector<double> r(3);
    double b = 0.0;
    double F = 0.0;
    for (int i = 0; i < DIM; ++i) {
        r[i] = atom1->pos[i] - atom2->pos[i];
        b += r[i] * r[i];
    }
    b = sqrt(b);
    F = (12*C12/pow(b, 14) - 6*C6/pow(b, 8));
    for (int i = 0; i < DIM; ++i) {
        r[i] /= b;
        atom1->force[i] += F * r[i];
        atom2->force[i] -= F * r[i];
    }
    return;
}

double exercises::System::bond_energy_x2(Bond *bond) {
    vector<double> r(3);
    double b = 0.0;
    for (int i = 0; i < DIM; ++i) {
        r[i] = bond->atom1->pos[i] - bond->atom2->pos[i];
        b += r[i] * r[i];
    }
    b = sqrt(b);
    return bond->k * pow((b - bond->b0), 2);
}

double exercises::System::lj_energy(Atom *atom1, Atom *atom2) {
    vector<double> r(3);
    double b = 0.0;
    for (int i = 0; i < DIM; ++i) {
        r[i] = atom1->pos[i] - atom2->pos[i];
        b += r[i] * r[i];
    }
    b = sqrt(b);
    return (C12/pow(b, 12) - C6/pow(b, 6));
}

void exercises::System::kick(double dt) {
    for (size_t i = 0; i < mols.size(); ++i) {
        // bond force
        for (Bond *bond : mols[i].bonds) {
            bond_force(bond);
        }
        // intermolecular LJ-force
        for (size_t j = i + 1; j < mols.size(); ++j) {
            for (Atom *atom1 : mols[i].atoms) {
                for (Atom *atom2 : mols[j].atoms) {
                    lj_force(atom1, atom2);
                }
            }
        }
        // propagate force into velocity and reset force
        for (Atom *atom : mols[i].atoms) {
            for (int i = 0; i < DIM; ++i) {
                atom->vel[i] += atom->force[i] * dt / atom->mass;
                atom->force[i] = 0.0;
            }
        }
    }
    return;
}

void exercises::System::drift(double dt) {
    for (Molecule &mol : mols) {
        for (Atom *atom : mol.atoms) {
            for (int i = 0; i < DIM; ++i) {
                atom->pos[i] += atom->vel[i] * dt;
            }
        }
    }
    return;
}

double exercises::System::kinetic_energy() {
    double E_kin = 0.0;
    double vel2;
    for (Molecule &mol : mols) {
        for (Atom *atom : mol.atoms) {
            vel2 = 0.0;
            for (int i = 0; i < DIM; ++i) {
                vel2 += atom->vel[i] * atom->vel[i];
            }
            E_kin += vel2 * atom->mass;
        }
    }
    return E_kin / 2.0;
}

double exercises::System::potential_energy() {
    double E_pot = 0.0;
    for (size_t i = 0; i < mols.size(); ++i) {
        // bond energy
        for (Bond *bond : mols[i].bonds) {
            E_pot += bond_energy_x2(bond);
        }
        E_pot /= 2.0;
        // intermolecular LJ-energy
        for (size_t j = i + 1; j < mols.size(); ++j) {
            for (Atom *atom1 : mols[i].atoms) {
                for (Atom *atom2 : mols[j].atoms) {
                    E_pot += lj_energy(atom1, atom2);
                }
            }
        }
    }
    return E_pot;
}

void exercises::System::print_atoms(int step, double timestep) {
    // better use a general output stream and/or use a buffer!
    // anyway stdout is buffered, isn't it? So that's fine for the moment.
    /*
     * Units:
     *   - time [ps]
     *   - distance / position [nm] --> print in Angstroms
     *   - energy [?] --> depend on the constants used
     *   - force [?] --> depend on the constants used
     */
    double energy = kinetic_energy() + potential_energy();
    step += 1;
    double time = step * timestep;
    cout.setf(ios::fixed, ios::floatfield);
    cout << "i = " << setw(8) << right << step << ", ";
    cout << "time = " << setw(8) << right << setprecision(3) << time << ", ";
    cout << "energy = " << setw(8) << right << setprecision(4) << energy << '\n';
    cout.setf(ios::scientific, ios::floatfield);
    for (Molecule &mol : mols) {
        for (Atom *atom : mol.atoms) {
            // print in Angstroms
            cout << setw(5) << left << "CH3";
            cout << setw(14) << setprecision(5) << left << 10*atom->pos[0];
            cout << setw(14) << setprecision(5) << left << 10*atom->pos[1];
            cout << setw(14) << setprecision(5) << left << 10*atom->pos[2] << '\n';
        }
    }
    return;
}

void exercises::System::md(int nsteps, double timestep) {
    double timestep_2 = timestep / 2.0;
    for (int i = 0; i < nsteps; ++i) {
        drift(timestep_2);
        kick(timestep);
        drift(timestep_2);
        print_atoms(i, timestep);
    }
    return;
}
